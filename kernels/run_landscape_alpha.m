function E = run_landscape_alpha
  sigma_range = -2:.25:2;
  s = load('images/PLANES', 'masks');
  
  %- form shape model
  [phi_basis phi_sigma phi_mean] = shape_model(s.masks);
  n = length(phi_sigma);
  
  %- prepare modeling function
  model_shape_fn = model_shape(phi_basis, phi_mean);
  
  %- form test image and reference model q
  [img h_fg h_bg] = test_image(s.masks);
  alpha = zeros([n 1]); % alpha=0 is mean shape
  q = model_shape_fn(img, [0 0], alpha);

  %- compute similarity across domain
  alpha = zeros([n 1]);
  rr = sigma_range*phi_sigma(1);
  cc = sigma_range*phi_sigma(2);
  for i = 1:length(rr)
    for j = 1:length(cc)
      alpha(1:2) = [rr(i) cc(j)];
      p = model_shape_fn(img, [0 0], alpha);
      E(i,j) = sum(sqrt(p.*q));
    end
  end
  
  %- display results
  clf;
  subplot(2,2,1); imagesc(img); axis image off;
  subplot(2,2,2); contourf(E, 20), axis ij square off;
  xx = 1:length(h_fg);
  subplot(2,2,3); plot(xx, h_fg, 'r', xx, h_bg, 'b');
  set(gca, 'XLim', [1 length(h_fg)], 'YTick', []);
  axis square; legend('fg', 'bg');
  subplot(2,2,4); plot(1:(256/length(q)):256,q,'r');
  set(gca, 'XLim', [1 length(h_fg)], 'YTick', []);
  axis square; legend('reference');
end




function [img h_fg h_bg] = test_image(masks)
  fg = struct('mu', 1/3*255, 'std', 20); % foreground
  bg = struct('mu', 2/3*255, 'std', 20); % background

  mu = zeros(size(masks{1}));
  for i = 1:length(masks)
    mu(:,:,2) = masks{i};
    mu = sum(mu, 3);
  end
  mu_mask = (mu / length(masks)) >= 0.5; % in at least half of masks?
  
  img =    bg.std*randn(size(masks{1})) + bg.mu;
  img_fg = fg.std*randn(size(masks{1})) + fg.mu;
  img(mu_mask) = img_fg(mu_mask);
  
  h_fg = imhist(uint8(img( mu_mask))); h_fg = h_fg / sum(h_fg);
  h_bg = imhist(uint8(img(~mu_mask))); h_bg = h_bg / sum(h_bg);
end


function [phi_basis phi_sigma phi_mean] = shape_model(masks)
  for i = 1:length(masks)
    x_train(:,i) = f(masks{i});
  end
  shape_fn = shape_lpca(x_train);
  [phi_basis phi_sigma phi_mean] = shape_fn.get_basis();
  
%   x = x_train(:,1);
%   phi_sigma'
%   alpha = phi_basis'*(x - mean(x_train,2))
%   x_ = shape_fn.preimage_alpha(alpha);
%   subplot(2,1,1); imagesc(f_(x));
%   subplot(2,1,2); imagesc(f_(x_));
%   keyboard

%   function m = f_(x)
%     m = reshape(x, size(masks{1}));
%   end

  function x = f(m)
    %x = single(m(:)); % mask
    x = reshape(mask2kernel(m), [numel(m) 1]); % kernel
  end
end
