function fn = model_tensor(kern)
  if ~exist('kern'), kern = @kern_Epanechnikov; end
  
  fn = @model;
  
  function T = model(img, center)
    %% calculate kernel
    win = (size(img) - 1)/2;
    [cc rr] = meshgrid((center(2)-win(2)):(center(2)+win(2)), ...
                       (center(1)-win(1)):(center(1)+win(1)));
    h = norm(max(size(img))/2);
    d2 = (cc.^2 + rr.^2)/h^2; % |(x-y)/h|^2
    K = kern(d2);
    
    %% compute feature vector for each point
    Z = feature_cov(img);

    %% compute weighted feature vectors
    z = Z(1,:);
    T = K(1)*z'*z;
    for i = 2:numel(img)
      z = Z(i,:);
      T = T + K(i)*z'*z;
    end
    T = T / sum(K(:)); % normalize kernel
  end
end


function k = kern_Epanechnikov(x)  % k' = constant
  d = 2; c_d = pi; % circle in 2D
  k = zeros(size(x));
  ind = find(x < 1);
  k(ind) = (d+2)*(1-x(ind))/(2*c_d);
end
