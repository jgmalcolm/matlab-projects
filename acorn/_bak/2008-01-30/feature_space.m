function fn = feature_space(win, b, m, epsilon)
  % win -- window dimensions [r c]
  % h -- scaling [r c]
  % b: I -> {1...m} pixel quantized in feature space
  % m -- number of quanta
  % epsilon -- Dirac delta epsilon

  d = 2; % circle
  c_d = pi;
  
  function k = kern(x) % Epanechnikov
    k = zeros(size(x));
    ind = find(x <= 1);
    k(ind) = (d+2)*(1-x(ind))/(2*c_d);
  end
  
  % pre-compute kernel window
  function [K C] = kern_mask(h)
    if ~exist('h','var'), h = 1; end
    scale = 1/norm(win/h)^2;
    [rr cc] = meshgrid(-win(2):win(2), -win(1):win(1));
    K = kern(scale*(rr.^2 + cc.^2)/h);
    C = 1/sum(K(:)); % normalization
  end
  
  [K1 C1] = kern_mask(1);

  function p = feature_vector(img, h)
    K = K1; C = C1; % default: no scaling
    if exist('h','var'), [K C] = kern_mask(h); end
    for i = 1:m
      p(i) = C * sum(sum(K.*dirac(b(img) - i, epsilon)));
    end
  end

  fn = @feature_vector;
end
