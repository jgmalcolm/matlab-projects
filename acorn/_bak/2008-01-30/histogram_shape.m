function fn = histogram_shape(shape, params)
  bins = params.bins;

  fn = @model;

  bin = @(x) round(x/(256/bins)); % which bin in feature space?
  dirac_lookup = dirac(0:bins+1, params.dirac_eps); % pre-compute
  k = size(shape.K_j, 3);
  function [p K K_j e_j d2_j U] = model(p_img, y, alpha)
    %- shift kernels
    K_j = zeros([size(p_img) k]);
    for j = 1:k
      K_j(:,:,j) = subpixel_shift(shape.K_j(:,:,j), y);
    end
    
    %- shape kernel directly from shape samples
    K = zeros(size(p_img));
    [d2_j e_j] = deal(zeros([1 k]));
    for j = 1:k
      d2_j(j) = (norm(alpha - shape.alpha_j(:,j))/shape.sigma)^2;
      e_j(j) = exp( -d2_j(j) );
      K = K + e_j(j)*K_j(:,:,j);
    end
    K = K / sum(K(:));
    
    %- sifting matrix [numel(p_img) bins]
    u = repmat(1:bins, [numel(p_img) 1]);
    B = repmat(bin(p_img(:)), [1 bins]);
    U = dirac_lookup(abs(B - u) + 1);
    
    %- form weighted histogram
    p = U' * K(:);
  end
end
