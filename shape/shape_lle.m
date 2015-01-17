function fn = shape_lle(X, method, k)
% SHAPE_LLE Shape analysis via local linear embedding

  fn.preimage = @preimage;
  
  n = size(X, 2);
  if ~exist('k'), k = ceil(n/3); end

  % kernels
  if exist('method')
    switch method
     case 'l', kern_init = @kern_linear;
     case 'k', kern_init = @kern_exp;
     case 'p', kern_init = @kern_poly;
    end
  else kern_init = @kern_linear; end
  kern = kern_init(X); % initialize

  %- precompute pairwise distances
  K = zeros(n);
  for i = 1:n
    for j = i:n
      [K(i,j) K(j,i)] = deal( kern(X(:,i), X(:,j)) ); % symmetric
    end
  end



  function x_hat = preimage(x)
    %- compute distance from each training sample
    for i = 1:n, k_x(i,1) = kern(x, X(:,i)); end
    k_xx = kern(x,x);
    %- gather neighborhood
    [d2 N] = sort(k_x, 'descend');
    N = N(2:(k+1));
    k_x = k_x(N);
    K_ = K(N,N);
    %- solve weights
    Q = zeros(k);
    for i = 1:k
      for j = 1:k
        [Q(i,j) Q(j,i)] = deal( k_xx - k_x(i) - k_x(j) + K_(i,j) );
      end
    end
    R = inv(Q);
    w = sum(R,2) / sum(R(:));
    
    %- reconstruct as weighted combination
    x_hat = X(:,N) * (w .* k_x) / dot(w, k_x);
  end
end









function [kern Dkern] = kern_exp(x_train)
  n = size(x_train, 2);

  %- determine variance in training data to set sigma
  dist = 0;  % method from Kwok's paper
  for i = 1:n
    for j = 1:n
      dist = dist + norm(x_train(:,i) - x_train(:,j))^2;
    end
  end
  sigma = dist / n / (n-1);

  kern = @(a, b) exp( -norm(a-b)^2 / sigma );
  Dkern = @(a,b,Da) -kern(a,b)/sigma * Da.*(a - b);
end



function [kern Dkern] = kern_linear(x_train)
  kern = @(a,b) dot(a, b);
  Dkern = @(a,b,Da) Da.*b;
end


function [kern Dkern] = kern_poly(x_train)
  c = 1; d = 3;
  kern = @(a,b) (c + dot(a,b))^d;
  Dkern = @(a,b,Da) 2*d*kern(a,b)/(c + dot(a,b)) * Da.*b;
end
