function fn = shape_kpca(x_train, norm2)
% SHAPE_KPCA Shape-based projection.
%
%  norm -- norm to use (default: L2 on difference)
%
%  >> fn = shape_kpca(x_train);
%  >> speed = fn.speed(x);
%  >> x_hat = fn.preimage(x);

  fn.speed = @speed;
  fn.preimage = @preimage;

  [d n] = size(x_train);

  k = min(7, n); % use at most seven eigenmodes

  % default norm is L2 distance.
  if ~exist('norm2','var'), norm2 = @(a,b) norm(a-b)^2; end
  
  % determine variance in training data
  dist = [];
  for i = 1:n
    for j = (i+1):n
      dist(end+1) = norm2(x_train(:,i), x_train(:,j));
    end
  end
  sigma = sqrt(mean(dist));
  kern = @(a, b)  exp( -norm2(a, b) / (2*sigma^2) );

  %%== (1) compute kernel K pairwise on training examples
  %%== (2) compute centered kernel K_tilde = H*K*H
  %                  where H = I - (1*1')/n and 1 = ones(1,n)
  %%== (3) SVD on K_tilde --> U*S*U'
  %%== (4) compute M = sum(u_i * u_i' / lambda_i)
  %%== (5) generate speed closure
  %%== (6) generate pre-image closure


  %%== (1) compute kernel K pairwise on training examples
  for i = 1:n
    for j = i:n
      K(i,j) = kern(x_train(:,i), x_train(:,j));
      K(j,i) = K(i,j); % symmetric
    end
  end
  
  %%== (2) compute centered kernel K_tilde = H*K*H
  H = eye(n) - 1/n*ones(n);
  K_tilde = H*K*H;

  %%== (3) SVD on K_tilde --> U*S*U'
  [U S U_T] = svd(K_tilde);
  lambda = diag(S);
  
  %%== (4) compute M = sum(u_i * u_i' / lambda_i)
  M = zeros(n);
  for i = 1:k
    M = M + U(:,i) * U(:,i)' / lambda(i);
  end



  %%== (5) generate closure
  function s = speed(x) % specific to exp(-norm(a-b)/sigma^2)
    %-- compute speed
    k = k_tilde(x);
    g = -2/n*ones(1,n) + 2*k'*M*K_tilde*H - 4*k'*M*H;
    s = zeros([d 1]);
    for i = 1:n
      x_i = x_train(:,i);
      s = s + g(i)*kern(x, x_i)*(x - x_i);
    end
    s = s / sigma^2;
    %-- homogenize with square distance (i.e. image term)
    s = 2*sigma^2*log((2 - s)/2);
  end
  




  %%== (6) pre-image closure
  function x_hat = preimage(x)
    %-- compute kernel of x against training set
    k_x = k_tilde(x);
    %-- compute d^2 distance of x against training set
    for i = 1:n
      k_xi = K(:,i);
      A = (k_x' + (2*K*ones(n,1)/n - 2*k_xi)'*H')*M*k_x;
      B = ones(1,n)*K*ones(n,1)/n^2;
      C = -2*ones(1,n)*k_xi/n;
      d2(i,1) = A + B + K(i,i) + C;
    end
    
    %-- beta projection coefficients
    beta = U(:,1:k)' * k_x ./ sqrt(lambda(1:k));
    
    %-- gamma weighting coefficients
    gamma = U(:,1:k) * beta;
    gamma(gamma < 0) = 0; % HACK: eliminate
    gamma = gamma + (1 - sum(gamma))/n; % normalize

    %-- form pre-image
    [num den] = deal(0);
    for i = 1:n % OPT: vectorize this
      coef = gamma(i) * (2 - d2(i))/2;
      num = num + coef * x_train(:,i);
      den = den + coef;
    end
    x_hat = num / den;
  end

  %%== (7) pre-image along mode of variation
  function x_hat = preimage_mode(mode, scale)
    error('unfinished');
    x_train_hat = U(:,mode)*U(:,mode)' * K_tilde / sqrt(lambda(mode));
    x_mean = mean(x_train_hat, 2);
    x_hat = x_mean + scale*std(x_train_hat,1,2);
  end



  function k = k_tilde(x) %%-- centered kernel computations
    for i = 1:n
      k(i,1) = kern(x, x_train(:,i));
    end
    k = H*(k - K*ones(n,1)/n); % center
  end

end
