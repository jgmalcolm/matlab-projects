function fn = filter_ukf(f_fn, h_fn, Q, R)
% FILTER_UKF Unscented Kalman filter for nonlinear systems.
%
% FN = FILTER_UKF(F, H, Q, R) Create an unscented Kalman filter:
%    F -- model transition function
%    H -- observation function
%    Q -- injected model noise
%    R -- injected observation noise
%
% Example: (1D point mass)
%  >> f_fn = @(x) [1 1; 0 1]*x;  % constant velocity
%  >> h_fn = @(x) [1 0]*x;       % position
%  >> x = [0 0]';
%  >> P = .1*eye(2);
%  >> est = filter_ukf(f_fn, h_fn, .1*eye(2), .1);
%  >> [x P] = est(x,P,4);  % x=[2.7 1.3]
%  >> [x P] = est(x,P,7);  % x=[6.4 2.5]
%  >> [x P] = est(x,P,10); % x=[9.8 3.0]
  
  fn = @filter;
  
  n = length(Q); % state dimension

  % Merwe's thesis, p52
  k = .01;
  w = [k; .5*ones(2*n,1)]/(n+k); w_ = diag(w);
  %w = [k .5]/(n+k);
  scale = sqrt(n + k);
  
  function [x P] = filter(x, P, z)
    %%-- (1) create sigma points
    X = sigma_points(x, P, scale);
    
    %%-- (2) predict state (also its mean and covariance)
    X = f_fn(X);
    x_hat = X * w;
    X_ = X - x_hat(:,ones(1,2*n+1)); % center
    P = X_ * w_ * X_' + Q; % partial update

    %%-- (3) predict observation (also its mean and covariance)
    Y = h_fn(X);
    Y_hat = Y * w;
    Y_ = Y - Y_hat(:,ones(1,2*n+1)); % center
    Pyy = Y_ * w_ * Y_' + R;

    %%-- (4) predict cross-correlation between state and observation
    Pxy = X_ * w_ * Y_';

    %%-- (5) Kalman gain K, estimate state/observation, compute cov.
    K = Pxy / Pyy;
    P = P - K*Pyy*K';
    x = x_hat + K*(z - Y_hat);
  end
end

function X = sigma_points(x, P, scale)
  M = scale * chol(P)';  % faster than matrix root since symmetric pos def
  X = x(:,ones(1,numel(x)));
  X = [x X+M X-M];
end
