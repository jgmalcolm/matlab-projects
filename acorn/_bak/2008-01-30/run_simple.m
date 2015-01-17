function run
  param.bins = 256;
  param.dirac_eps = 10;
  param.mu = 0.3;
  param.max_iterations = 10;

  paths; clf; colormap gray;
  %s = load('images/SYNTH_constellation');
  %param.win = s.win;
  s = load('images/3BB', 'imgs');
  param.win = 5;

  %- reference density q
  [hist_fn K dK] = histogram(param);
  img  = extract(s.imgs{1}, s.win, round(s.init));
  y = s.init;
  q = hist_fn(double(img));

  for i = 1:length(s.imgs)
    y = acorn(y, q, s.imgs{i}, s.win, dK, hist_fn, param, i);
  end
end





function y_ = acorn(y, q, img, win, dK, hist_fn, param, id)
  bin = @(x) round(x/(256/param.bins)); % which bin in feature space?
  dirac_lookup = dirac(0:param.bins+1, param.dirac_eps); % pre-compute

  [cc rr] = meshgrid(-param.win(2):param.win(2), -param.win(1):param.win(1));
  cc = dK(:).*cc(:); rr = dK(:).*rr(:);

  for t = 1:param.max_iterations
    %- candidate patch
    p_img = double(extract(img, param.win, round(y)));
    [p U] = hist_fn(p_img);
    
    %- 1. weights from similarity metric
    d_dp = sqrt(q./(p+eps));
    w = U * d_dp;
    
    %- 2. factor out y
    dy = sum([w.*rr w.*cc])'/(sum(w)+eps);
    y_ = y + dy;

    %- 3. display
    subplot(2,2,1);
    imagesc(img); axis image off;
    title(['frame=' int2str(id) '     iter=' int2str(t)]);
    hold on; plot(y_(2), y_(1), 'r*'); hold off;

    subplot(2,2,2);
    imagesc(p_img); axis image off;

    subplot(2,2,3);
    imagesc(reshape(w, size(p_img))); axis image off;
    title('w(i)');
    
    subplot(2,2,4);
    xx = 1:length(q);
    plot(xx,p,'r', xx,q,'b');
    
    drawnow;

    if norm(y_ - y) < param.mu, break, end
    y = y_;
  end
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
