function h = warping_speed
  % Dambreville et al. "Shape-based approach to robust image segmentation using kPCA" CVPR 2006

  % parameters
  param.eps = 4;

  % required
  h.init_iteration = @init_iteration;
  h.move_in = @(p) 0; % dummy
  h.move_out = @(p) 0; % dummy
  % extensions
  h.init = @init;
  
  s = [];  % dummy for lexical scoping


  function init(masks)
    s.fn = shape_pca(cell2mat(map(@(x) double(x(:)), masks)), 'k');
  end

  
  E = [];
  function S = init_iteration(phi, C)
    HP = heaviside(phi(:), param.eps);
    dP = dirac(phi(:), param.eps);
    
    Dd2 = s.fn.gradient_d2(HP, dP);
    
    S = Dd2(C) - 0*2*kappa(phi, C);

    sp(2,1,1); imagesc(zeros(size(phi)), [0 1]); axis image off
    hold on; contour(phi, [0 0], 'w', 'LineWidth', 2); hold off;
    drawnow
    
    E(end+1) = s.fn.d2_x_Px(HP);
    sp(2,1,2); plot(E);
  end
end





function H = heaviside(x, eps)
  H = zeros(size(x));
  H(eps < x) = 1;
  ind = find(abs(x) <= eps);
  x = x(ind);
  H(ind) = (1 + x/eps + sin(pi*x/eps)/pi)/2;
end

function d = dirac(x, eps)
  d = zeros(size(x));
  ind = find(abs(x) <= eps);
  x = x(ind);
  d(ind) = (1 + cos(pi*x/eps))/2/eps;
end
