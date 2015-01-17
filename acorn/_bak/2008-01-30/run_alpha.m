function run_alpha(ref)
  paths;
  s = load('images/SYNTH4_PLANE');
  if ~exist('ref'), ref = 1; end

  %- form shape model
  [phi_basis phi_sigma phi_mean k shape_fn] = shape_model(s.masks);

  %- prepare reference image and model q
  [img alpha_ref dy] = test_image(s.masks{ref}, shape_fn);
  model_fn = model_shape(phi_basis, phi_mean);
  q = model_fn(img, [0 0], alpha_ref);

  weight_fn = weight_histogram(q);
  
  clf;

  y = s.win + dy; % center
  alpha = zeros(size(alpha_ref)); % mean
  [y alpha t] = acorn_shape(img, s.win, model_fn, weight_fn, y, alpha, ...
                              phi_basis, phi_sigma);
  
  fprintf('start: %f\nfinal: %f\n', norm(alpha_ref), norm(alpha_ref - alpha));
end



function [img alpha dy] = test_image(mask, shape_fn)
  fg = struct('mu', 255/3,   'std', 10);
  bg = struct('mu', 2*255/3, 'std', 10);
  
  y = weighted_centroid(mask);
  dy = round(y) - y;

  x = reshape(mask2kernel(mask), [numel(mask) 1]);
  alpha = shape_fn.get_alpha(x);
  
  img =    fg.std*randn(size(mask)) + fg.mu;
  img_bg = bg.std*randn(size(mask)) + bg.mu;
  img(mask==1) = img_bg(mask==1);
end



function [phi_basis phi_sigma phi_mean k shape_fn] = shape_model(masks)
  for i = 1:length(masks)
    x_train(:,i) = f(masks{i});
  end
  shape_fn = shape_lpca(x_train); % linear PCA
  [phi_basis phi_sigma phi_mean] = shape_fn.get_basis();
  k = length(phi_sigma);

  function x = f(m)
%     m = imdilate(m, [0 1 0; 1 1 1; 0 1 0]);
%     m = imdilate(m, [0 1 0; 1 1 1; 0 1 0]);
    x = reshape(mask2kernel(m), [numel(m) 1]); % kernel
  end
end
