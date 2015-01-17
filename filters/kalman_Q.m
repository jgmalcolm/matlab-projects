function err = kalman_Q(A, B, C, Q, R, x0, tr)
%KALMAN_Q Calculate L2-norm error based on Q used.  See filter_kalman.m
%
% FN = FILTER_KALMAN(A, B, C, Q, R, P0, XO, TR) Create a Kalman filter for the
% linear system:
%    dx/dt = A x + B u
%        y = C x
% Q describes system noise, R describes measurement noise, initial state is
% X0, and TR is the ground truth trackpoints.
%
% Example: (1D point mass)
%  >> err = kalman_Q([1 1;0 1], [0 0]', [1 0], ...
%                    zeros(2), .1, [0 0]', ...
%                    [[0 0]' [1 1]' [2 2]']);
  
  fn = filter_kalman(A, B, C, Q, R, x0);
  
  for i=1:length(tr)
    y = fn(tr(:,i), 0);
    err(i) = norm(y - tr(:,i));
  end
