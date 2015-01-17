function run
  paths;
  param.iterations = 80;
  param.write_frame = @() 1; %write_frame();
  %param.appearance = @appearance_mean;
  %param.appearance = @appearance_croppedmean;
  param.appearance = @appearance_original;
  %param.appearance = @appearance_nobackground;

  L = param;
  L.dt_x = 2e-6;
  K = param;
  K.dt_x = 2e-6;
  
  clf; colormap gray;
  %s = load('images/WALKING'); s.masks = s.masks_L_; s.init  = s.init_L; % left person
  %s = load('images/SHARK'); s.masks = s.masks_full;
  s = load('images/SOCCER'); s.masks = s.masks_full;

%   train_masks = {s.masks_R_{1:3:end} s.masks_L_{1:3:end}};
%   train_imgs  = {s.imgs{1:3:end} s.imgs{1:3:end}};
%   s.init = weighted_centroid(train_masks{1});
  
  [L.win K.win] = deal(s.win);
  train_masks = {s.masks{1:4:end}};
  train_imgs  = {s.imgs{1:4:end}};
  [L.shape_fn L.ss_plot] = linear_shape_model(train_masks, train_imgs, L);
  [K.shape_fn K.ss_plot] = kernel_shape_model(train_masks, train_imgs, K);

  [L.x K.x] = deal(s.init(:));
  for i = 1:length(s.masks);
    L.x = ktrack(L.x, s.imgs{i}, L.shape_fn, L, L.ss_plot, 1);
    K.x = ktrack(K.x, s.imgs{i}, K.shape_fn, K, K.ss_plot, 2);
    disp(' ');
  end
end


function x = ktrack(x, img_full, shape_fn, param, ss_plot, fig)
  for t = 1:param.iterations
    dx = x - round(x);
    %- reference
    img = double(extract(img_full, param.win+1, round(x)));
    img = subpixel_shift(img, -dx); % register image to atlas
    img = extract(img, param.win, param.win+2);

    %- preimage and gradients of
    [img_hat E_x alpha] = shape_fn(img, dx);
    %- solve for x
    x = x - param.dt_x*E_x;
    
    %- display
    sp(2,1,fig);
    subplot(1,4,1); imagesc(img, [0 255]); axis image off;
    hold on; plot(param.win(2), param.win(1), 'r+'); hold off;
    title(['image (' int2str(t) ')']);
    subplot(1,4,2); imagesc(img_hat, [0 255]); axis image off;
    title('pre-image');
    img(img_hat == 0) = 0;
    subplot(1,4,3); imagesc((img - img_hat).^2); axis image off;
    title('squared difference');
    subplot(1,4,4); ss_plot(alpha); axis square;
    title('first two modes of shape space');
    drawnow; param.write_frame();
  end
  fprintf('E=%.0f\n', sum((img(:) - img_hat(:)).^2));
end


function [preimage_fn ss_plot] = linear_shape_model(masks, imgs, param)
  preimage_fn = @preimage;
  ss_plot = @plot_shape_space;
  
  % form training data
  win = param.win;
  for i = 1:length(masks)
    c = round(weighted_centroid(masks{i}));
    mask = extract(masks{i}, win, c);
    img = extract(imgs{i}, win, c);
    x_train(:,i) = param.appearance(mask, img);
  end
  sz = size(img); % size of windowed image
  
  % form shape basis via linear PCA
  shape_fn = shape_lpca(x_train);
  
  for i = 1:length(masks)
    alpha_j(:,i) = shape_fn.get_alpha(x_train(:,i));
  end
  function plot_shape_space(a)
    plot(alpha_j(1,:), alpha_j(2,:), 'k.', a(1), a(2), 'ro');
  end
  
  % sample image on uniform grid (while we may center kernel off grid)
  [cc rr] = meshgrid(-param.win(2):param.win(2), -param.win(1):param.win(1));
  %- grid coordinates
  xx = [rr(:)'; cc(:)']; xx(3,:) = 1;
  
  % computes preimage and its gradients wrt x and alpha
  function [img_hat E_x alpha] = preimage(img, dx)
    %- preimage
    alpha = shape_fn.get_alpha(img(:));
    img_hat = reshape(shape_fn.preimage_alpha(alpha), sz);
    
    %- grad_x preimage
    x_x(:,1) = reshape(Dy(img_hat), [numel(img_hat) 1]); % note: y is 1st
    x_x(:,2) = reshape(Dx(img_hat), [numel(img_hat) 1]);
    dX = x_x'; dX(3,:) = 0;
    [Ty Tx] = T_partials; % transform partials
    P(1,:) = sum(dX .* (Ty * xx)); % y translation
    P(2,:) = sum(dX .* (Tx * xx)); % x translation
    %img(img_hat == 0) = 0;
    E_x = P * (img(:) - img_hat(:));
  end
end



function [preimage_fn ss_plot] = kernel_shape_model(masks, imgs, param)
  preimage_fn = @preimage;
  ss_plot = @plot_shape_space;

  % form training data
  win = param.win;
  for i = 1:length(masks)
    c = round(weighted_centroid(masks{i}));
    mask = extract(masks{i}, win, c);
    img = extract(imgs{i}, win, c);
    x_train(:,i) = param.appearance(mask, img);
  end
  sz = size(img); % size of windowed image
  
  % form shape basis via kernel PCA
  shape_fn = shape_pca(x_train);
  
  [U lambda K] = shape_fn.get_basis();
  
  for i = 1:length(masks)
    beta_j(:,i) = shape_fn.get_beta(x_train(:,i));
  end
  function plot_shape_space(b)
    plot(beta_j(1,:), beta_j(2,:), 'k.', b(1), b(2), 'ro');
  end
  
  % sample image on uniform grid (while we may center kernel off grid)
  [cc rr] = meshgrid(-param.win(2):param.win(2), -param.win(1):param.win(1));
  %- grid coordinates
  xx = [rr(:)'; cc(:)']; xx(3,:) = 1;
  
  % computes preimage and its gradients wrt x and alpha
  function [img_hat E_x beta] = preimage(img, dx)
    %- preimage
    [beta k_x] = shape_fn.get_beta(img(:));
    img_hat = reshape(shape_fn.preimage_beta(beta, k_x), sz);
    
    %- grad_x (distance in input space)
    x_x(:,1) = reshape(Dy(img_hat), [numel(img_hat) 1]); % note: y is 1st
    x_x(:,2) = reshape(Dx(img_hat), [numel(img_hat) 1]);
    dX = x_x'; dX(3,:) = 0;
    [Ty Tx] = T_partials; % transform partials
    P(1,:) = sum(dX .* (Ty * xx)); % y translation
    P(2,:) = sum(dX .* (Tx * xx)); % x translation
    %img(img_hat == 0) = 0;
    E_x = P * (img(:) - img_hat(:));
  end
end



function [Ty Tx] = T_partials
  Tx = [0 0 0; 0 0 1; 0 0 0];
  Ty = [0 0 1; 0 0 0; 0 0 0];
end


function x = appearance_mean(m, img)
  u = mean(img(m));
  v = mean(img(~m));
  x = u*m + v*~m;
  x = x(:);
end
function x = appearance_original(m, img)
  x = double(img(:));
end

function x = appearance_nobackground(m, img)
  x = zeros([numel(img) 1]);
  x(m) = img(m);
end

function x = appearance_croppedmean(m, img)
  u = mean(img(m));
  m_out = m;
  m_out = imdilate(m_out, [0 1 0; 1 1 1; 0 1 0]);
  m_out = imdilate(m_out, [0 1 0; 1 1 1; 0 1 0]);
  m_out = imdilate(m_out, [0 1 0; 1 1 1; 0 1 0]);
  m_out = imdilate(m_out, [0 1 0; 1 1 1; 0 1 0]);
  m_out = m_out & ~m;
  v = mean(img(m_out));
  x = u*m + v*~m;
  x = x(:);
end
