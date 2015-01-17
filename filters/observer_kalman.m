function [est_fn pred_fn] = observer_kalman(A, B, C, D, Q, R, x)
% OBSERVER_KALMAN Kalman observer for linear systems.
%
% FN = OBSERVER_KALMAN(A, B, C, D, Q, R, XO) Create a Kalman observer for the
% linear system:
%    dx/dt = A x + B u
%        y = C x + D u
% Q describes system noise, R describes measurement noise, and initial state
% is X0.
%
% Example: (1D point mass)
%  >> fn = observer_kalman([1 1;0 1], [0 0]', [1 0], [0 0]', ...
%                          zeros(2), .1, [0 0]');
%  >> y = fn(4); % 3.8095
%  >> y = fn(7, 0); % 6.8421 -- input optional
  
  est_fn = @estimate;
  pred_fn = @predict;

  p_hat = ones(numel(x)); % TODO: arbitrary?
  
  last_z = [];
  last_u = [];

  function x_tild = predict(z, u)
    if ~exist('z','var'), z = z_last; u = u_last; end
    x_tild = A*x + B*u;
  end

  function y = estimate(z, u)
    if ~exist('u','var'), u = zeros([size(B,2) size(x,2)]); end
    last_z = z; last_u = u;
    x_tild = predict(z, u);
    p_tild = A*p_hat*A' + Q;
    K = (p_tild * C') / (C*p_tild*C' + R);
    x = x_tild + K*(z - C*x_tild);
    p_hat = (eye(size(K*C)) - K*C)*p_tild;
    y = C*x + D*u;
  end
  
end
