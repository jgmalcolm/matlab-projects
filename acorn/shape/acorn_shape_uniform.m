function [y alpha t] = acorn_shape(img, win, model_fn, weight_fn, y, alpha, phi_basis, phi_sigma)
  N_max = 40; % max iterations
  dt_y = 10;
  dt_alpha = phi_sigma / 100;

  % where sample image on uniform grid (while we may center kernel off grid)
  [cc rr] = meshgrid(-win(2):win(2), -win(1):win(1));
  %- grid coordinates
  xx = [cc(:)'; rr(:)']; xx(3,:) = 1;

  for t = 1:N_max
    %- candidate patch
    p_img = single(extract(img, win, round(y)));
    dy = y - round(y);
    [p K] = model_fn(p_img, dy, alpha);
    
    %- weights from similarity metric
    w = weight_fn(p_img, p);

    %%-- y
    % kernel derivatives
    Kx = Dx(K); Ky = Dy(K);
    dK = [Ky(:)'; Kx(:)']; dK(3,:) = 0;
    [Ty Tx] = T_partials; % transform partials
    % weighted gradient
    drho_y(1) = sum(dK .* (Ty * xx)) * w;
    drho_y(2) = sum(dK .* (Tx * xx)) * w;
    y = y - dt_y*drho_y;
    
    %%-- alpha
    phi_basis_ = shift_basis(dy, phi_basis, 2*win+1);
    drho_alpha = w' * phi_basis_;
    alpha = alpha + dt_alpha.*drho_alpha';
  end
end




function phi_basis = shift_basis(dx, phi_basis, sz)
  %- shift basis
  n = size(phi_basis,2);
  phi_basis = reshape(phi_basis, [sz n]);
  for i = 1:n
    phi_basis(:,:,i) = subpixel_shift(phi_basis(:,:,i), dx);
  end
  phi_basis = reshape(phi_basis, [prod(sz) n]);
end



function [Ty Tx] = T_partials
  Tx = [0 0 0; 0 0 1; 0 0 0];
  Ty = [0 0 1; 0 0 0; 0 0 0];
end
