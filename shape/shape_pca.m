function fn = shape_pca(x_train, method)
% SHAPE_PCA Shape analysis via kernel PCA (default: linear)
%
%  >> fn = shape_pca(x_train);
%  >> alpha = [0 .2 .2 0 -.5];
%  >> x_hat = fn.preimage_alpha(alpha);

  fn.get_basis = @get_basis;
  fn.get_beta = @get_beta;
  fn.preimage = @preimage;
  fn.preimage_beta = @preimage_beta;
  fn.gradient_d2 = @gradient_d2;
  fn.d2_x_Px = @d2_x_Px;

  % form basis from eigenvectors comprising top xx% of training data
  param.basis_threshold = 0.99;

  n = size(x_train, 2);
  
  % kernels (default: exponential)
  if exist('method') && method(1) == 'l', kern_init = @kern_linear;
  elseif exist('method') && method(1) == 'k', kern_init = @kern_exp;
  elseif exist('method') && method(1) == 'p', kern_init = @kern_poly;
  else kern_init = @kern_exp; end

  %%== (1) compute kernel K pairwise on training examples
  %%== (2) compute centered kernel K_tilde
  %%== (3) SVD on K_tilde --> U*S*U'
  %%== (4) determine basis size based on eigenvalue contribution
  %%== (5) form reduced orthonormal basis
  %%== (6) generate closures: x_hat, alpha, modes, speed


  %%== (1) compute kernel K pairwise on training examples
  [kern Dkern] = kern_init(x_train); % initialize kernel
  for i = 1:n
    for j = i:n
      [K(i,j) K(j,i)] = deal( kern(x_train(:,i), x_train(:,j)) ); % symmetric
    end
  end

  %%== (2) compute centered kernel K_tilde = H*K*H
  I = ones(n,1);
  H = eye(n) - 1/n*ones(n); % centering matrix
  K_tilde = H*K*H;

  %%== (3) SVD on K_tilde --> U*S*U'
  [U S] = svd(K_tilde);
  sigma = sqrt(diag(S));
  
  %%== (4) determine basis size based on eigenvalue contribution
  if ~exist('k')
    % grab top xx% eigenvalues
    sigma_ = cumsum(sigma/sum(sigma));
    k = find(sigma_ >= param.basis_threshold, 1, 'first');
  end
  
  %%== (5) form truncated orthonormal basis and covariance M
  sigma = sigma(1:k);
  U = U(:,1:k);
  M = zeros(n);
  for i = 1:k % reduced size
    M = M + U(:,i) * U(:,i)' / sigma(i);
  end


  %%== (6) generate closures
  
  function [U_ sigma_ K_] = get_basis
    U_ = U;
    sigma_ = sigma;
    K_ = K;
  end

  function [beta k_x] = get_beta(x)  % projection coefficients
    for i = 1:n, k_x(i,1) = kern(x, x_train(:,i)); end
    beta = U' * H * (k_x - K*I/n);
  end

  %- pre-image from back-projection (note: assumes EXP kernel?)
  function x_hat = preimage(x)
    [beta k_x] = get_beta(x); % projection coefficients
    x_hat = preimage_beta(beta, k_x);
  end
  
  function x_hat = preimage_beta(beta, k_x)
    %-- gamma weighting coefficients
    gamma = U * beta; % unnormalize U
    gamma(gamma < 0) = 0; % HACK: stabilize
    gamma = gamma + (1 - sum(gamma))/n; % recenter

    %-- compute d^2 distance of x against training set
    for i = 1:n
      k_xi = K(:,i);
      k_xi_tilde = H*(K(:,i) - K*I/n);
      k_x_tilde  = H*(k_x    - K*I/n);
      d2(i) = (k_x_tilde - 2*k_xi_tilde)'*M*k_x_tilde ...
              + I'*K*I/n^2 + K(i,i) - 2*I'*k_xi/n;
    end
    
    %-- form pre-image as weighted combination
    [num den] = deal(0);
    for i = 1:n % OPT: vectorize this
      coef = gamma(i) * (2 - d2(i))/2;
      num = num + coef * x_train(:,i);
      den = den + coef;
    end
    x_hat = num / den;
  end

  %- pre-image along mode of variation
  function x_hat = preimage_mode(mode, scale)
    beta = zeros(n,1);
    beta(mode) = scale * sqrt(sigma(mode));
    x_hat = preimage_beta(beta);
  end

  %- gradient: D_x d^2(x,Px)
  function Dx = gradient_d2(x, dx)
    dK_x = zeros(numel(x), n);
    for i = 1:n
      dK_x(:,i) = Dkern(x, x_train(:,i), dx);
    end

    for i = 1:n, k_x(i,1) = kern(x, x_train(:,i)); end
    k_x_tilde = H*(k_x - K*I/n);
    
    Dx = - 2*dK_x*(k_x_tilde'*M*H + I'/n)' + Dkern(x,x,dx);
  end
  
  function d = d2_x_Px(x)
    for i = 1:n, k_x(i,1) = kern(x, x_train(:,i)); end
    k_x_tilde = H*(k_x - K*I/n);
    d = -k_x_tilde' * M * k_x_tilde + I'*K*I/n^2 + kern(x,x) -2*I'*k_x/n;
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

  kern = @(a,b) exp( -norm(a-b)^2 / sigma );
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
