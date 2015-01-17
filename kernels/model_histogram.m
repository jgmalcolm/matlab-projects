function fn = model_histogram(sz)
  bins = 64;
  epsilon = 1;

  fn = @model;

  bin = @(x) ceil(x/(256/bins)); % which bin in feature space?
  dirac_lookup = dirac(0:bins-1, epsilon); % precompute Dirac (max_dist=bins-1)
  
  %- pre-compute kernel (assume centered at origin)
  win = (sz - 1)/2;
  [cc rr] = meshgrid(-win(2):win(2), -win(1):win(1));
  h = norm(win);
  K = kern(cc(:)/h, rr(:)/h);
  K = K / sum(K);

  function [p dK] = model(p_img, offset)
    %- compute sifting matrix [numel(p_img) bins]
    B = repmat(bin(p_img(:)), [1 bins]);
    u = repmat(1:bins, [numel(p_img) 1]); % OPT: pre-compute
    U = dirac_lookup(abs(B - u) + 1);
    %U = (B == u); % HACK: unsmoothed
    
    %- kernel derivatives
    [dK{1} dK{2}] = deal(1);

    %- form histogram
    p = U' * K;
  end
end



function k = kern(xx, yy)  % Epanechnikov => constant k'
  d2 = xx.^2 + yy.^2; % |x-y|^2
  d = 2; c_d = pi; % circle in 2D
  d2(1 < d2) = 1;
  k = (d+2)*(1-d2)/(2*c_d);
end
