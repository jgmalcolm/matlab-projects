function fn  = shape_lpca(x_train)
% SHAPE_LPCA Shape analysis via linear PCA.
%
%  >> fn = shape_lpca(x_train);
%  >> alpha = [0 .2 .2 0 -.5]; % assume top 5 eigenvalues used
%  >> x_hat = fn.preimage_alpha(alpha);

  %fn.get_basis_cnt = @get_basis_cnt;
  fn.get_basis = @get_basis;
  fn.get_alpha = @get_alpha;
  fn.preimage_mode = @preimage_mode;
  fn.preimage_alpha = @preimage_alpha;
  fn.preimage = @preimage;

  % form basis from eigenvectors comprising top 90% of training data
  param.basis_threshold = 0.99;

  [d n] = size(x_train);
  
  %%== (1) compute mean shape (x_mean)
  %%== (2) compute mean offset maps (M_i)
  %%== (3) SVD:  M'*M/n = U*Q*U' (using kernel trick)
  %%== (4) determine basis size based on eigenvalue contribution
  %%== (5) form orthonormal basis
  %%== (6) generate closures: x, alpha

  
  %%== (1) compute mean shape (x_mean)
  x_mean = mean(x_train, 2);

  %%== (2) compute mean offset maps
  M  = x_train - x_mean(:, ones(1,n));

  %%== (3) SVD:  M'*M/n = U*Q*U'
  [U S V] = svd(M' * M);  % unnormalized U
  sigma = sqrt(diag(S));
  
  %%== (4) determine basis size based on eigenvalue contribution
  sigma_ = cumsum(sigma/sum(sigma));
  k = find(sigma_ >= param.basis_threshold, 1, 'first');
  %%== (5) form orthonormal basis: U_i = M'*U_i
  U = M * U(:,1:k);
  sigma = sigma(1:k);
  for i = 1:k
    U(:,i) = U(:,i) * sigma(i);
  end


  %%== (6) generate closures

  function x_hat = preimage_mode(mode, scale)
    x_hat = U(:,mode)*scale + x_mean;
  end

  %- pre-image from backprojection
  function x_hat = preimage(x)
    alpha = get_alpha(x); % projection cooefficients
    x_hat = preimage_alpha(alpha);
  end

  %- pre-image from alpha
  function x_hat = preimage_alpha(alpha)
    x_hat = U*alpha + x_mean;
  end
  
  %- alpha coefficients
  function a = get_alpha(x)
    a = U'*(x - x_mean);
  end


  %- get basis size
  function k_ = get_basis_cnt
    k_ = k;
  end
  
  %- get basis
  function [U_ sigma_ x_mean_] = get_basis
    U_ = U;
    sigma_ = sigma;
    x_mean_ = x_mean;
  end

end
