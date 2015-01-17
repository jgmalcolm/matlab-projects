function [y alpha] = acorn_shape(y, alpha, img, win, hist_fn, params)

  % sample image on uniform grid (while we may center kernel off grid)
  [cc rr] = meshgrid(-win(2):win(2), -win(1):win(1));
  %- grid coordinates
  xx = [rr(:)'; cc(:)']; xx(3,:) = 1;

  bin = @(x) round(x/(256/bins)); % which bin in feature space?
  dirac_lookup = dirac(0:bins+1, params.dirac_eps); % pre-compute

  for t = 1:params.iterations
    %- candidate patch
    p_img = single(extract(img, win, round(y)));
    dy = y - round(y);
    [p K K_set x_set alpha_set] = hist_fn(p_img, dy, alpha);
    
    %- weights from similarity metric
    w = weight_fn(p_img, p);
    v = shape_weights(K_set, x_set, w, params);

    %%-- y   (gradient descent)
    % kernel derivatives
    Kx = Dx(K); Ky = Dy(K);
    dK = [Ky(:)'; Kx(:)']; dK(3,:) = 0;
    [Ty Tx] = T_partials; % transform partials
    % weighted gradient
    drho_y(1) = sum(dK .* (Ty * xx)) * w;
    drho_y(2) = sum(dK .* (Tx * xx)) * w;
    y = y - params.dt_dy*drho_y;
    
    %%-- alpha   (fixed point iteration)
    alpha = alpha_set * v / (sum(v(:)) + eps);
    
    %- display
    subplot(2,2,1);
    imagesc(overlay(uint8(p_img), K>0.00025)); axis image off;
    title(['iter=' int2str(t)]);

    subplot(2,2,2);
    a = 1; b = 2;
    is_counted = 0 <= x_set & x_set < 1;
    plot(alpha_set(a,is_counted), alpha_set(b,is_counted), 'b.', ...
         alpha_set(a,~is_counted), alpha_set(b,~is_counted), 'k.', ...
         alpha(a), alpha(b), 'r*');
    axis square
    title(sprintf('shape space: %dx%d', a, b));

    subplot(2,2,3);
    %imagesc(reshape(w, size(p_img))); axis image off;
    w_ = reshape(w, size(p_img));
    imagesc(w_); axis image off; title('w(i)');
    
    subplot(2,2,4);
    x = 1:length(x_set);
    plot(x, v/sum(v), 'b', x, x_set, 'r', x, [1 length(x_set), 'k:');
    legend('v(j)', '|\alpha-\alpha_j|', 'Location', 'BestOutside');
    
    keyboard
    drawnow;
  end
end





function [Ty Tx] = T_partials
  Tx = [0 0 0; 0 0 1; 0 0 0];
  Ty = [0 0 1; 0 0 0; 0 0 0];
end


function v = shape_weights(K_, x, w, params)
  % prepare derivatives (HACK: approximate)
  dKe = zeros(size(x)); % default: zero for |a-a_j|/sigma > 1
  dKe(0 <= x & x < 1) = 1; % otherwise: assume constant
  % prepare K_j set for easy multiplication
  sz = size(K_);
  K_ = reshape(K_, [prod(sz(1:2)) sz(3)]);
  % done
  v = dKe' .* (K_' * w);
  return

  sz = size(K_);
  K_ = reshape(K_, [prod(sz(1:2)) sz(3)]);
  v = K_' * w;
  [v_ ind] = sort(v');
  v(ind(1:end-params.k)) = 0;
  return
end
