function speculate(imgs, pts)
  param.bins = 128;
  param.dirac_eps = 4;
  param.win = [12 12];

  paths; clf; colormap gray;
  %s = load('images/rob', 'imgs', 'pts');
  %s = load('images/marvin_and_bottle', 'imgs', 'pts');
  n = min(numel(pts), 3);

  %- reference density q
  hist_fn = histogram(param);
  xx = 1:param.bins;

  for t = 20:5:length(imgs)
    subplot(4,n,1:3*n);
    imagesc(imgs{t}); axis image off; title(t);

    for i = 1:n
      y = pts{i}(t,:);
      p_img = double(extract(imgs{t}, param.win, round(y)));
      h = hist_fn(p_img);

      %subplot(4,n,1:3*n); hold on; plot(y(2), y(1), 'r.'); hold off;
      subplot(4,n,1:3*n); hold on; plot_box(y, param.win, 'r'); hold off;
      subplot(4,n,3*n+i); plot(xx, h); set(gca, 'YTickLabel', []);
    end
    drawnow;
  end
end




function plot_box(y, w, c)
  corners = [y + [-1 -1].*w; ...
             y + [-1 +1].*w; ...
             y + [+1 +1].*w; ...
             y + [+1 -1].*w; ...
             y + [-1 -1].*w];
  plot(corners(:,2), corners(:,1), c);
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
