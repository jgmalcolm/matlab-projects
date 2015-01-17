function run
  param.bins = 64;
  param.dirac_eps = 3;
  param.mu = 0.05;
  param.max_iterations = 15;
  param.win = [9 9];

  param.scale = 10;
  param.Q = 0.1*eye(2);
  param.R = 0.1*eye(2);

  paths; clf; colormap gray;
  s = load('images/marvin_and_bottle_2', 'imgs', 'init');

%   imagesc(s.imgs{1}); axis image off;
%   [cc rr] = getpts();
  %y = s.init(:,3) - [4 0]';
  y = s.init(:,2);

  %- reference density q
  [hist_fn K dK] = histogram(param);
  img  = extract(s.imgs{1}, param.win, round(y));
  q = hist_fn(double(img));

  %- set up UKF
  % prediction: identity
  f_fn = @(x, u) x;
  % measurement: mean-shift
  h_fn = @(y, img) acorn(y, q, img, dK, hist_fn, param);
  filter = filter_ukf(f_fn, h_fn, param.Q, param.R, y, param.scale);

  for t = 1:20 %length(s.imgs)
%     subplot(3,2,1:4);
    imagesc(s.imgs{t}); axis image off; title(t); hold on;
%     y = acorn(y, q, s.imgs{t}, dK, hist_fn, param);
    y = filter(s.imgs{t});
    plot(y(2), y(1), 'g.');
%     subplot(3,2,1:4);
    hold off;
    drawnow;
  end
end





function y_ = acorn(y, q, img, dK, hist_fn, param)
  bin = @(x) round(x/(256/param.bins)); % which bin in feature space?
  dirac_lookup = dirac(0:param.bins+1, param.dirac_eps); % pre-compute
  win = param.win;
  [cc rr] = meshgrid(-win(2):win(2), -win(1):win(1));
  cc = dK(:).*cc(:); rr = dK(:).*rr(:);

%   subplot(3,2,1:4);
  plot(y(2), y(1), 'bo');
  
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
