function run_tensor
  %s = load('images/ISRAEL_movie14');
  s = load('images/SYNTH2');

  % initial target
  target = single(extract(s.imgs{1}, s.win, round(s.init)));
  model_fn = model_tensor();
  q = model_fn(target, s.init - round(s.init));
  
  weight_fn = weight_tensor(q);
  
  figure(1); clf; colormap gray;
  subplot(2,2,1); imagesc(target); axis image; title('target');
  subplot(2,2,2); imagesc(q); axis image; title('target tensor');
  
  y = s.init;
  tic;
  for i = 1
    img = s.imgs{i};

    subplot(2,2,3); imagesc(img); axis image
    hold on; plot(y(2), y(1), 'ro'); hold off;

    [y t] = acorn(y, img, s.win, model_fn, weight_fn);

    hold on; plot(y(2), y(1), 'r.'); hold off;
    drawnow;
  end
  fprintf('%.2f Hz\n', length(s.imgs)/toc);

end
