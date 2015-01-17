function run
  param.bins = 256;
  param.dirac_eps = 1;
  param.mu = 0.05;
  param.max_iterations = 15;
  param.win = 7*[1 1];

  param.is_filtered = false; % turn filter on/off
  param.scale = 5;
  param.Q = .5*eye(4);
  param.R = .1*eye(4);

  paths; clf; colormap gray;
  s = load('images/marvin_and_bottle_2', 'imgs', 'init');

%   x = s.init(:,1);
%   x = s.init(:,2);
  x = s.init(:,3) - [4 0]';
  x(3:4) = 0;

  %- reference density q
  [hist_fn K dK] = histogram(param);
  img  = extract(s.imgs{1}, param.win, round(x(1:2)));
  q = hist_fn(double(img));
  
%   subplot(2,1,1); imagesc(img);
%   subplot(2,1,2); plot(q);
%   return

  %- set up UKF
  %- model: idenity
  %f_fn = @(x, u) x;
  %h_fn = @(x, img) acorn(x, q, img, dK, hist_fn, param);
  %- model: constant velocity
  f_fn = @(x, u) [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1] * x;
  % measurement: mean-shift
  function [y e] = h_fn(x, img)
    [y e] = acorn(x(1:2), q, img, dK, hist_fn, param);
    y(3:4,1) = x(3:4);
  end
  filter = filter_ukf(f_fn, @h_fn, param.Q, param.R, x, param.scale);
  
  xx = 1:numel(q);

  for t = 1:200 %length(s.imgs)
    subplot(4,1,1:3);
    imagesc(s.imgs{t}); axis image off; title(t); hold on;
    if param.is_filtered,  x = filter(s.imgs{t});
    else                   x(1:2) = acorn(x(1:2), q, s.imgs{t}, dK, hist_fn, param);
    end
    X(:,t) = x(1:2);
    plot(X(2,:), X(1,:), 'y', 'LineWidth', 2);
    hold off;

    subplot(4,1,4);
    img  = extract(s.imgs{t}, param.win, round(x(1:2)));
    p = hist_fn(double(img));
    plot(xx,q,'r',  xx,p,'b');

    drawnow;
  end
end





function [y_ e] = acorn(y, q, img, dK, hist_fn, param)
  win = param.win;
  [cc rr] = meshgrid(-win(2):win(2), -win(1):win(1));
  cc = dK(:).*cc(:); rr = dK(:).*rr(:);

%   subplot(3,2,1:4);
  plot(y(2), y(1), 'bo');
  box = [y' + [-1 -1].*win; ...
         y' + [-1  1].*win; ...
         y' + [ 1  1].*win; ...
         y' + [ 1 -1].*win; ...
         y' + [-1 -1].*win];
  plot(box(:,2), box(:,1), 'g');
  
  xx = 1:numel(q);

  for t = 1:param.max_iterations
    %- candidate patch
    p_img = double(extract(img, win, round(y)));
    [p U] = hist_fn(p_img);
    
%     subplot(3,2,5); plot(xx',p,'r', xx',q,'b');
    
    %- 1. weights from similarity metric
    d_dp = sqrt(q./(p+eps));
    w = U * d_dp;
    
%     subplot(3,2,6); imagesc(reshape(w, size(p_img))); axis image off

    %- 2. factor out y
    dy = sum([w.*rr w.*cc])'/(sum(w)+eps);
    y_ = y + dy;

    %- 3. test convergence
    if norm(y_ - y) < param.mu, break, end
    y = y_;
  end
%   subplot(3,2,1:4);
  plot(y_(2), y_(1), 'r.');

  e = sum(sqrt(p .* q));
end








function [fn K dK] = histogram(params)
  bins = params.bins;

  fn = @model;

  bin = @(x) round(x/(256/bins)); % which bin in feature space?
  dirac_lookup = dirac(0:bins+1, params.dirac_eps); % pre-compute

  %- pre-compute kernel (assume centered at origin)
  win = params.win;
  [cc rr] = meshgrid(-win(2):win(2), -win(1):win(1));
  h = min(win);
  [K dK] = kern(cc(:)/h, rr(:)/h);
  K = K / sum(K);
  
  function [p U] = model(p_img)
    %- sifting matrix [numel(p_img) bins]
    u = repmat(1:bins, [numel(p_img) 1]);
    B = repmat(bin(p_img(:)), [1 bins]);
    U = dirac_lookup(abs(B - u) + 1);
    
    %- form weighted histogram
    p = U' * K(:);
  end
end


function [K dK] = kern(xx, yy)  % Epanechnikov
  d2 = xx.^2 + yy.^2; % |x-y|^2
  d = 2; c_d = pi; % circle in 2D
  d2(1 < d2) = 1; % limit support

  K = (d+2)*(1-d2)/(2*c_d);

  dK = (d+2)/(2*c_d) * ones(size(K));
  dK(1 <= d2) = 0;
end



function d = dirac(x, eps)
%   d = zeros(size(x));
%   ind = find(abs(x) < eps);
%   x = x(ind);
%   d(ind) = (1 + cos(pi*x/eps))/2/eps;
  d = eps./(eps^2 + x.^2)/pi;  % infinite support
end





function fn = filter_ukf(f_fn, h_fn, Q, R, x, k)
% FILTER_UKF Unscented Kalman filter for nonlinear systems.
%
% FN = FILTER_UKF(F, H, Q, R, P0, X0) Create an unscented Kalman filter using
% model F, observation H, system noise Q, measurement noise R, and initial
% covariance P0 and state X0.
%
% Example: (1D point mass)
%  >> f_fn = @(x,u) [1 1; 0 1]*x;
%  >> h_fn = @(x,u) [1 0]*x;
%  >> fn = filter_ukf(f_fn, h_fn, zeros(2), .1, [0 0]');
%  >> y = fn(4, 0); % 3.8
%  >> y = fn(7, 0); % 6.8
  
  fn = @filter;

  if ~exist('k'), k = 1; end % sigma-point scaling parameter
  n = length(x);
  P = eye(n);
  w = [k/(n + k)  1/2/(n + k)];

  function x_ = filter(u)
    %%-- (1) create sigma points
    X = sigma_points(x, P, k);

    %%-- (2) predict state (also its mean and covariance)
    for i = 1:2*n+1
      X(:,i) = f_fn(X(:,i), u); %-- states
    end
    x_hat = w_mean(X, w); %%-- mean
    Xm = X - x_hat(:,ones(1,2*n+1)); % mean-shifted
    P = w_cov(Xm, Xm, w) + Q; %%-- covariance

    %%-- (3) predict observation (also its mean and covariance)
    for i = 1:2*n+1
      [Y(:,i) e(i)] = h_fn(X(:,i), u); %-- observations
    end
    Y_hat = w_mean(Y, w);
    Ym = Y - Y_hat(:,ones(1,2*n+1)); % mean-shifted
    Pyy = w_cov(Ym, Ym, w) + R;

    %%-- (4) predict cross-correlation between state and observation
    Pxy = w_cov(Xm, Ym, w);

    %%-- (5) Kalman gain K, estimate state/observation, compute cov.
    K = Pxy / Pyy;
    P = P - K*Pyy*K';
    [max_e max_id] = max(e);
    z = Y(:,max_id);
    plot(z(2), z(1), 'g*');
    x = x_hat + K*(z - Y_hat);
    x_ = x;
  end
end

function X = sigma_points(x, P, k)
  n = length(x);
  X(:,1) = x;
  M = sqrtm((n+k)*P);
  for i = 1:n
    X(:,1+i)   = x + M(:,i);
    X(:,1+i+n) = x - M(:,i);
  end
%   X = x(:,ones([1 2*n+1])) + [0*x M -M]; % vectorized
end

function c = w_cov(X, Y, w) %-- weighted covariance
  n = size(X,1);
  c = w(1) * X(:,1) * Y(:,1)';
  for i = 1:n
    c = c + w(2)*(X(:,1+i)*Y(:,1+i)' + X(:,1+i+n)*Y(:,1+i+n)');
  end
end

function m = w_mean(X, w) %-- weighted mean
  n = size(X,1);
  m = X * w([1 2*ones(1,2*n)])';
end
