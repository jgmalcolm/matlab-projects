function C = model_cov(img, center, kern)
  if ~exist('kern'), kern = @kern_Epanechnikov; end
  N = numel(img);
  
  %% determine kernel weights
  rows = size(img,1); cols = size(img,2);
  [cc rr] = meshgrid((1:cols) - (cols+1)/2 - center(2), ...
                     (1:rows) - (rows+1)/2 - center(1));
  h = norm(max(size(img))/2);
  d2 = (cc.^2 + rr.^2)/h^2; % |(x-y)/h|^2
  K = kern(d2);
  
  %% compute feature vector for each point
  Z = feature_cov(img);

  %% compute weighted covariances
  z = Z(1,:);
  C = K(1)*z'*z;
  for i = 2:N
    z = Z(i,:);
    C = C + K(i)*z'*z;
  end
  C = C / sum(K(:)); % normalize

  %% compute weighted mean feature vector
  mu = K(1)*Z(1,:);
  for i = 2:N
    mu = mu + K(i)*Z(i,:);
  end
  mu = mu / sum(K(:));
  
  %% subtract off mean
  C = C - mu' * mu;
end


function k = kern_Epanechnikov(x)  % k' = constant
  d = 2; c_d = pi; % circle in 2D
  k = zeros(size(x));
  ind = find(x < 1);
  k(ind) = (d+2)*(1-x(ind))/(2*c_d);
end
