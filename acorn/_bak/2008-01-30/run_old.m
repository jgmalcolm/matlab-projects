function run
  
  data = load('images/ISRAEL_movie14');
  data.imgs = data.imgs_color; % use color
  
  epsilon = 8;
  
  % initial target
  target = extract(data.imgs{1}, data.win, data.init);
  q = fn(single(target));
  
  figure(1); clf; colormap gray;
  subplot(2,2,1); imagesc(target); title('target');
  subplot(2,2,2); plot(q); title('target feature profile');
  
  tr = [];
  y = data.init;
  N = size(data.imgs, 3);
  N = 8;
  tic;
  for i = 1:N
    img = data.imgs(:,:,i);

    subplot(2,2,3); imagesc(img); hold on;
    plot(y(2), y(1), 'y*');
    [y t] = acorn(y, img, data.win, b, m, epsilon, q, fn);
    plot(y(2), y(1), 'r*');

    subplot(2,2,4);
    img_win = extract(img, 1*data.win, round(y));
    L = landscape(img_win, data.win, q, fn);
    mesh(L);
    pt = repmat(data.win, [10 1]);
    hold on; plot3(pt(:,2), pt(:,1), .1:.1:1, 'r.'); hold off;

    drawnow;
    fprintf('[pause]...\n'); pause;
    tr(:,end+1) = [y t]';
  end
  T = toc;
  fprintf('%.2f Hz\n', N/T);

end




function fn = model_histogram(sz)
  bins = 64;
  epsilon = 1;

  fn = @model;

  bin = @(x) ceil(x/(256/bins)); % which bin in feature space?
  dirac_lookup = dirac(0:bins-1, epsilon); % precompute Dirac (max_dist=bins-1)
  
  %- pre-compute kernel (assume centered at origin)
  win = (sz - 1)/2;
  [cc rr] = meshgrid(-win(2):win(2), -win(1):win(1));
  h = norm(win);
  K = kern(cc(:)/h, rr(:)/h);
  K = K / sum(K);

  function h = model(p_img, offset)
    %- compute sifting matrix [numel(p_img) bins]
    B = repmat(bin(p_img(:)), [1 bins]);
    u = repmat(1:bins, [numel(p_img) 1]); % OPT: pre-compute
    U = dirac_lookup(abs(B - u) + 1);
    %U = (B == u); % HACK: unsmoothed
    
    %- form histogram
    h = U' * K;
  end
end



function k = kern(xx, yy)  % Epanechnikov => constant k'
  d2 = xx.^2 + yy.^2; % |x-y|^2
  d = 2; c_d = pi; % circle in 2D
  d2(1 < d2) = 1;
  k = (d+2)*(1-d2)/(2*c_d);
end
