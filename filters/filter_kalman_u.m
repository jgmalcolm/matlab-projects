function fn = filter_kalman(A, B, C, D, Q, R)
% FILTER_KALMAN Classic Kalman filter for linear systems.
%
% FN = FILTER_KALMAN(A, B, C, D, Q, R) Create a Kalman filter for the linear
% system:
%    dx/dt = A x + B u
%        y = C x + D u
% Q describes system noise, R describes measurement noise, and initial state
% is X0.
%
% Example: (1D point mass)
%  >> fn = filter_kalman([1 1;0 1], [0 0]', [1 0], [0 0]', ...
%                        zeros(2), .1, [0 0]');
%  >> y = fn(4); % 3.8095
%  >> y = fn(7, 0); % 6.8421 -- input optional
  
  p_hat = ones(length(x)); % TODO: arbitrary?
  
  function y = filter(z, u)
    if ~exist('u'), u = zeros([size(B,2) size(x,2)]); end
    x_tild = A*x + B*u;
    p_tild = A*p_hat*A' + Q;
    K = (p_tild * C') / (C*p_tild*C' + R);
    x = x_tild + K*(z - C*x_tild);
    p_hat = (eye(size(K*C)) - K*C)*p_tild;
    y = C*x + D*u;
  end
  
  fn = @filter; % return closure reference
end
