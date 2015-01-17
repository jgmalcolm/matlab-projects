%============================= dfilter_kalman ============================
%
% DFILTER_KALMAN Classic Kalman filter for discrete linear systems.
%
% dkfilt = DFILTER_KALMAN(A, B, C, Q, R, XO, P0) 
%
% Create a Kalman filter for the linear system:
%
%    dx/dt = A x + B u
%        y = C x
%
% Q describes system noise, R describes measurement noise, and 
% the initial state and covariance are is X0 and P0, respectively.
%
% Example: (1D point mass, unit time, no forcing)
%
%  >> dkfilt = dfilter_kalman([1 1;0 1], [], [1 0], ...
%                        zeros(2), .1, [0 1]', [0.5 0;0 0.5]);
%
%  ... process loop ...
%
%  >> xpred = dkfilt.predict();
%
%  ...
%
%  >> [y, x] = dkfilter.correct(ymeasured); 
%
%  ...
%
%
%  Copyright 2007.
%
%============================= dfilter_kalman ============================

%
%  Name:	dfilter_kalman.m
%
%  Author:	Jimi Malcolm,		(Created)
%		Patricio A. Vela.	(Modified)
%
%  Created:	XX/XX/2006
%  Modified:	07/13/2007
%
%============================= dfilter_kalman ============================
function dkfilt = dfilter_kalman(A, B, C, Q, R, x, P)
  
  dkfilt.predict = @kf_predict;
  dkfilt.correct = @kf_correct;

  x_hat = x;
  P_hat = P;

  function [xpred, Ppred] = kf_predict
    x_hat = A*x_hat;
    P_hat = A*P_hat*A' + Q;

    if (nargout == 2)
      xpred = x_hat;
      Ppred = P_hat;
    elseif (nargout == 1)
      xpred = x_hat;
    end
  end

  function [y, x] = kf_correct(z, u)
    if (nargin == 2), x_tild = x_tild + B*u; end
    K = (P_hat * C')*inv((C*P_hat*C' + R));
    x_hat = x_hat + K*(z - C*x_hat);
    P_hat = (eye(size(P_hat)) - K*C)*P_hat;
    y = C*x_hat;
    if (nargout == 2)
      x = x_hat;
    end
  end
  
end
