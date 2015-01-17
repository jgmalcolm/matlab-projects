% TODO: ensure shapes are properly aligned for [0 0]
function fn = model_shape(K_set)
  bins = 64;
  dirac_epsilon = 1;

  fn = @model;

  bin = @(x) round(x/(256/bins)); % which bin in feature space?
  dirac_lookup = dirac(0:bins+1, dirac_epsilon); % precompute Dirac
  
  function [p K K_ Ke alpha_set] = model(p_img, offset, alpha)
    %- shift kernels
    for j = 1:size(K_set.K,3)
      K_(:,:,j) = subpixel_shift(K_set.K(:,:,j), offset);
    end
    
    %- shape kernel directly from shape samples
    K = zeros(size(K_set.K(:,:,1)));
    for j = 1:size(K_,3)
      Ke(j) = kern(norm(alpha - K_set.alpha(:,j))/K_set.sigma);
      K = K + Ke(j)*K_(:,:,j);
    end
    K = K / sum(Ke(:));
    alpha_set = K_set.alpha;
    
    keyboard

    subplot(2,1,2);
    h = plot(K_set.alpha(1,:), K_set.alpha(2,:), 'k.', ...
             alpha(1), alpha(2), 'r*');
    drawnow;

    %- sifting matrix [numel(p_img) bins]
    u = repmat(1:bins, [numel(p_img) 1]);
    B = repmat(bin(p_img(:)), [1 bins]);
    U = dirac_lookup(abs(B - u) + 1);
    
    %- form weighted histogram
    p = U' * K(:);
  end
end


function k = kern(x) % Epanechnikov (k' => constant)
  % TODO: do we need the d,c_d scalars?
  d = 2; c_d = pi; % circle in 2D
  ind = find(0 <= x & x < 1);
  k = zeros(size(x));
  k(ind) = (d+2)*(1-x(ind))/(2*c_d);
end
