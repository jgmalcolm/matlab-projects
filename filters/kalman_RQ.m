function [R Q] = kalman_RQ(f, g, x, u, z)
% KALMAN_RQ Maximizes joint likelihood of filter to determine R and Q.
%
% Naive approach from "Discriminative Training of Kalman Filters" P Abbeel, A
% Coates, M Montemerlo, A Ng, S Thurn.
  
  T = size(x,2);
  
  R = zeros(size(z,1));
  for i=1:T
    d = z(:,i) - g(x(:,i));
    R = R + d*d';
  end
  R = R / (T+1);

  Q = zeros(size(x,1));
  for i=2:T
    d = x(:,i) - f(x(:,i-1), u(:,i));
    Q = Q + d*d';
  end
  Q = Q / T;
