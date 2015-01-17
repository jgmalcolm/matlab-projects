function run_histogram
  
  s = load('images/SYNTH2');
  %s = load('images/SYNTH3_PLANE');
  
  % initial target
  target = single(extract(s.imgs{1}, s.win, round(s.init)));
  model_fn = model_histogram(size(target));
  q = model_fn(target, s.init - round(s.init));
  weight_fn = weight_histogram(q);
  
  % kernel
  dK_fn = dK_Epanechnikov(size(target));
  
  figure(1); clf; colormap gray;
  subplot(2,2,1); imagesc(target); title('target');
  xx = 1:256/length(q):256;
  subplot(2,2,2); plot(xx,q); title('target density profile');
  set(gca, 'XLim', [1 256]);
  
  y = s.init;
  tic;
  for i = 1:length(s.imgs)
    img = s.imgs{i};

    subplot(2,2,3); imagesc(img);
    hold on; plot(y(2), y(1), 'ro'); hold off;

    [y t] = acorn(img, s.win, model_fn, weight_fn, dK_fn, y);

    hold on; plot(y(2), y(1), 'r.'); hold off;
    drawnow;
  end
  fprintf('%.2f Hz\n', length(s.imgs)/toc);

end


function dK_fn = dK_Epanechnikov(sz)
  dK_fn = @(y0, alpha) 1;
end
