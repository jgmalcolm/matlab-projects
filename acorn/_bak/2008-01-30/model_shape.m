function fn = model_shape(K_set, params)
  bins = params.bins;

  fn = @model;

  bin = @(x) round(x/(256/bins)); % which bin in feature space?
  dirac_lookup = dirac(0:bins+1, params.dirac_eps); % pre-compute
  
  function [p K K_ x alpha_set] = model(p_img, offset, alpha)
    %- shift kernels
    for j = 1:size(K_set.K,3)
      K_(:,:,j) = subpixel_shift(K_set.K(:,:,j), offset);
    end
    
    %- shape kernel directly from shape samples
    K = zeros(size(K_set.K(:,:,1)));
    for j = 1:size(K_,3)
      x(j) = norm((alpha - K_set.alpha(:,j))*K_set.sigma);
      K = K + kern(x(j))*K_(:,:,j);
    end
    K = K / sum(K(:));
    %K = K / size(K_,3);
    alpha_set = K_set.alpha;
    
    %- sifting matrix [numel(p_img) bins]
    u = repmat(1:bins, [numel(p_img) 1]);
    B = repmat(bin(p_img(:)), [1 bins]);
    %U = dirac_lookup(abs(B - u) + 1);
    U = (B == u);
    
    %- form weighted histogram
    p = U' * K(:);
  end
end


function k = kern(x) % Epanechnikov (k' => constant)
  d = 2; c_d = pi; % assume 2D circle
  ind = find(0 <= x & x < 1);
  k = zeros(size(x));
  k(ind) = (d+2)*(1-x(ind))/(2*c_d);
end
