function run_fixed
  paths; clf; colormap gray;
  s = load('images/SYNTH4_PLANE');

  %- compute shape subspace and samples within
  K_set = shape_model(s.masks);

  %- prepare reference image and model q
  model_fn = model_shape(K_set);
  img  = extract(s.imgs{1}, s.win, round(s.init));
  y = s.init;
  q = model_fn(single(img), y-round(y), K_set.alpha(:,1));
  weight_fn = weight_histogram(q);
  
  %alpha = zeros(size(K_set.alpha(:,1))); % mean shape
  alpha = K_set.alpha(:,1);

  img = single(s.imgs{1});
  [y_ alpha_] = acorn_shape(img, s.win, model_fn, weight_fn, y, alpha);

  subplot(2,1,1); imagesc(img, [0 255]); axis image off;
  hold on; plot(y(2), y(1), 'r*'); hold off;

  %print('-dpng', sprintf('/tmp/out/frame_%03d.png', i));
  drawnow; fprintf('\n');

end


function K_set = shape_model(masks)
  % form training data
  for i = 1:length(masks)
    [x_train(:,i) K(:,:,i)] = f(masks{i});
  end
  
  % form shape basis
  shape_fn = shape_lpca(x_train); % linear PCA
  %shape_fn = shape_pca(x_train); % kernel PCA
  
  
  % determine alpha of each training sample
  for i = 1:length(masks)
    alpha(:,i) = shape_fn.get_alpha(x_train(:,i));
  end
  
  % determine max pairwise distance between alpha -> sigma
  dist = [];
  for i = 1:length(masks)
    for j = i+1:length(masks)
      dist(end+1) = norm(alpha(:,i) - alpha(:,j));
    end
  end
  sigma = mean(dist);
  
  K_set = struct('K', K, 'alpha', alpha, 'sigma', sigma);

  function [x kern] = f(m)
    %m = imdilate(m, [0 1 0; 1 1 1; 0 1 0]);
    kern = mask2kernel(m);
    x = reshape(kern, [numel(m) 1]); % kernel
  end
end
