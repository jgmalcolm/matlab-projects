function fn = histogram(params)
  bins = params.bins;

  fn = @model;

  bin = @(x) round(x/(256/bins)); % which bin in feature space?
  dirac_lookup = dirac(0:bins+1, params.dirac_eps); % pre-compute

  %- pre-compute kernel (assume centered at origin)
  win = params.win;
  [cc rr] = meshgrid(-win(2):win(2), -win(1):win(1));
  h = norm(win);
  K = kern(cc(:)/h, rr(:)/h);
  K = K / sum(K);

  function [p U] = model(p_img, y)
    %- shift kernel
    K_ = subpixel_shift(K, y);
    
    %- sifting matrix [numel(p_img) bins]
    u = repmat(1:bins, [numel(p_img) 1]);
    B = repmat(bin(p_img(:)), [1 bins]);
    U = dirac_lookup(abs(B - u) + 1);
    
    %- form weighted histogram
    p = U' * K_(:);
  end
end


function K = kern(xx, yy)  % Epanechnikov
  d2 = xx.^2 + yy.^2; % |x-y|^2
  d = 2; c_d = pi; % circle in 2D
  d2(1 < d2) = 1; % limit support
  K = (d+2)*(1-d2)/(2*c_d);
end
