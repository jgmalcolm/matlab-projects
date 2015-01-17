function run_shape
  param.bins = 128;
  param.dirac_eps = 5;
  param.alpha0 = 4;
  param.iterations = 40;
  param.sigma_scale = 0.5;
  param.dt_dy = 0.5;
  param.dt_dalpha = 0.1; % SHARK
  %param.dt_dy = 20;
  %param.dt_dalpha = 0.5; % SHARK
  %param.dt_dalpha = 0.1; % SOCCER
  %param.dt_dalpha = 0.01; % SYNTH5_shapes

  paths; clf; colormap gray;
  %s = load('images/SYNTH3_PLANE');
  %s = load('images/SYNTH5_shapes');
  s = load('images/SHARK');
  %s = load('images/SOCCER');
  %s = load('images/WALKING'); s.init  = s.init_L;
  %s.masks = {s.masks{1:3:end}}; % reduce training set

  %- compute shape subspace and samples within
  shape = form_shape_model(s.masks, param);

  %- reference density q
  hist_fn = histogram_shape(shape, param);
  img  = extract(s.imgs{1}, s.win, round(s.init));
  y = s.init;
  alpha = shape.alpha_j(:,1);
  q = hist_fn(single(img), y-round(y), alpha);
  
  alpha = shape.alpha_j(:,param.alpha0);
  for i = 1 %length(s.imgs)
    [y alpha] = acorn(y, alpha, q, s.imgs{i}, s.win, hist_fn, shape, param, i);
  end
end





function m = form_shape_model(masks, param)
  % form training data
  for i = 1:length(masks)
    [x_j(:,i) K_j(:,:,i)] = f(masks{i});
  end
  
  % form shape basis
  shape_fn = shape_lpca(x_j); % linear PCA
  %shape_fn = shape_pca(x_j); % kernel PCA
  
  % determine alpha of each training sample
  for j = 1:length(masks)
    alpha(:,j) = shape_fn.get_alpha(x_j(:,j));
  end
  
  % determine max pairwise distance between alpha -> sigma
  dist = [];
  for j = 1:length(masks)
    for k = j+1:length(masks)
      dist(end+1) = norm(alpha(:,j) - alpha(:,k));
    end
  end
  sigma = mean(dist)*param.sigma_scale;

  m = struct('K_j', K_j, 'alpha_j', alpha, 'sigma', sigma);

  function [x K_j] = f(m)
    m = imdilate(m, [0 1 0; 1 1 1; 0 1 0]);
    %m = imdilate(m, [0 1 0; 1 1 1; 0 1 0]);
    K_j = mask2kernel(m);
    x = reshape(K_j, [numel(m) 1]); % shape kernel
  end
end






function [y alpha] = acorn(y, alpha, q, img, win, hist_fn, shape, param, id)
  % sample image on uniform grid (while we may center kernel off grid)
  [cc rr] = meshgrid(-win(2):win(2), -win(1):win(1));
  %- grid coordinates
  xx = [rr(:)'; cc(:)']; xx(3,:) = 1;

  bin = @(x) round(x/(256/param.bins)); % which bin in feature space?
  dirac_lookup = dirac(0:param.bins+1, param.dirac_eps); % pre-compute

  for t = 1:param.iterations
    %- candidate patch
    p_img = double(extract(img, win, round(y)));
    dy = y - round(y);
    [p K K_j e_j d2_j U] = hist_fn(p_img, dy, alpha);
    
    %- 1. weights from similarity metric
%     d_dp = sqrt(q./(p+eps));
%     w = U * d_dp;
    w = load('images/SYNTH5_shapes', 'masks_full');
    w = load('images/SHARK', 'masks_full');
    w = double(extract(w.masks_full{id}, win, round(y)));
    w = w(:);
    
    %- 2. shape weights from kernels
    K_j = reshape(K_j, [numel(w) size(K_j,3)]);
    v = e_j' .* (K_j' * w);
    
    %- 3. solve for y
    Kx = Dx(K); Ky = Dy(K);    % kernel derivatives
    dK = [Ky(:)'; Kx(:)']; dK(3,:) = 0;
    Ty = [0 0 1; 0 0 0; 0 0 0]; % transform partial derivatives
    Tx = [0 0 0; 0 0 1; 0 0 0];
    drho_y(1) = sum(dK .* (Ty * xx)) * w; % weighted gradient
    drho_y(2) = sum(dK .* (Tx * xx)) * w;
    y = y - param.dt_dy*drho_y;
    
    %- 4. solve for alpha
    drho_alpha = zeros([numel(alpha) 1]);
    for j = 1:numel(v)
      drho_alpha = drho_alpha + (alpha - shape.alpha_j(:,j))*v(j);
    end
    alpha = alpha - param.dt_dalpha*drho_alpha;
    
    %- 5. display
    subplot(2,2,1);
    imagesc(overlay(uint8(p_img), K>0.00025)); axis image off;
    title(['frame=' int2str(id) '     iter=' int2str(t)]);

    subplot(2,2,2);
    a = 1; b = 2;
    plot(shape.alpha_j(a,:), shape.alpha_j(b,:), 'k.', ...
         shape.alpha_j(a,id), shape.alpha_j(b,id), 'bx', ...
         alpha(a), alpha(b), 'ro');
    axis square
    title(sprintf('shape space: %dx%d', a, b));

    subplot(2,2,3);
    imagesc(reshape(w, size(p_img))); colorbar; axis image off;
    title('w(i)');
    
    subplot(2,2,4);
    x = 1:length(e_j);
    plot(x, v/sum(v), 'b', x, e_j, 'r');
    
    drawnow;
  end
end
