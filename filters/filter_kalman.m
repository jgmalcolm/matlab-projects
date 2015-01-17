function fn = filter_kalman(F, H, Q, R)
% FILTER_KALMAN Classic Kalman filter for linear systems.
%
% FN = FILTER_KALMAN(F, H, Q, R) Create a Kalman filter for the linear system:
%    dx/dt = F x
%        y = H x
%    Q -- injected model noise
%    R -- injected observation noise
%
% Example: (1D point mass)
%  >> est = filter_kalman([1 1;0 1], [1 0], .1*eye(2), .1);
%  >> x = [0 0]';
%  >> P = .1*eye(2);
%  >> [x P] = est(x,P,  4); % x=[3.1 1.0]
%  >> [x P] = est(x,P,  7); % x=[6.4 2.2]
%  >> [x P] = est(x,P, 10); % x=[9.7 2.8]
  
  fn = @filter;

  function [x P] = filter(x, P, z)
    x_ = F*x;
    P_ = F*P*F' + Q;
    K = (P_ * H') / (H*P_*H' + R);
    x = x_ + K*(z - H*x_);
    P = P_ - K*H*P_;
  end
  
end
