function fn  = shape_lpca(x_train)
% SHAPE_LPCA Shape-based projection via linear PCA.
%
%  >> fn = shape_pca(x_train);
%  >> speed = fn.speed(x);
%  >> x_hat = fn.preimage(x);

  fn.speed = @speed;
  fn.preimage = @preimage;
  fn.preimage_mode = @preimage_mode;

  [d n] = size(x_train);

  k = min(7, n);

  %%== (1) compute mean shape (x_mean)
  %%== (2) compute mean offset maps (M_i)
  %%== (3) SVD:  M'*M/n = U*Q*U' (using kernel trick)
  %%== (4) orthonormal basis: J_i = M'*U_i / sqrt(q_i)
  %%== (5) generate speed closure
  %%== (6) generate pre-image closure
  %%== (6) generate pre-image closure along mode of variation

  
  %%== (1) compute mean shape (x_mean)
  x_mean = mean(x_train, 2);

  %%== (2) compute mean offset maps (x_off)
  M  = x_train - x_mean(:, ones(1,n));

  %%== (3) SVD:  M'*M/n = U*Q*U' (using transpose trick)
  [U_ S V] = svd(M' * M);  % unnormalized U
  lambda = diag(S);

  %%== (4) orthonormal basis: U_i = M'*U_i / sqrt(q_i)
  lambda = lambda(1:k);
  U = M * U_;
  U = U(:,1:k) ./ sqrt(lambda(:, ones(1,d))');
  

  %%== (5) generate closure
  function s = speed(x)
    x_hat = preimage(x);
    s = 2*(x - x_hat);
  end

  %%== (6) pre-image closure
  function x_hat = preimage(x)
    alpha = U'*(x(:) - x_mean); % projection coordinates
    x_hat = U*alpha + x_mean;
  end
  
  %%== (7) pre-image along mode of variation
  function x_hat = preimage_mode(mode, scale)
    x_hat = U(:,mode)*scale*sqrt(lambda(mode)) + x_mean;
  end
end
