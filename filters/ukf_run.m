function ukf_run
n=3;      %number of state
q=0.1;    %std of process 
r=0.1;    %std of measurement
Q=q^2*eye(n); % covariance of process
R=r^2;        % covariance of measurement  
f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
h=@(x)x(1)+x(2);                     % measurement equation

s=[0;0;1];                                % initial state
x=s+q*randn(3,1); %initial state          % initial state with noise
P = eye(n);                               % initial state covraiance
N=200;                                     % total dynamic steps
xV = zeros(n,N);          % Cao estmate
xV_= zeros(n,N);          % my UKF estmate
sV = zeros(n,N);          % actual
zV = zeros(1,N);

my_ukf = filter_ukf(f, h, @my_z, Q, R, x);

for k=1:N
  z = h(s) + r*randn;                     % measurments
  sV(:,k)= s;                             % save actual state
  zV(:,k)= z;                             % save measurment
  [x, P] = ukf(f,x,P,h,z,Q,R);            % ekf 
  xV(:,k) = x;                            % save estimate
  xV_(:,k) = my_ukf();
  s = f(s) + q*randn(3,1);                % update process 
end
for k=1:3                                 % plot results
  subplot(3,1,k)
  plot(1:N, sV(k,:), 'y-', ...
       1:N, xV(k,:),  'm--', ...
       1:N, xV_(k,:)+.002, 'r--');
end

xV = xV(:) / norm(xV(:));
xV_= xV_(:) / norm(xV_(:));
sV = sV(:) / norm(sV(:));
cao_err = sum((xV - sV).^2)
 my_err = sum((xV_ - sV).^2)

function z = my_z(x)
  z = zV(:,k);
end
end
