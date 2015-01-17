function E = landscape(img, win, alpha, model_fn, dist_fn)
  target = single(img);

  y0 = 0;
  alpha = 0;
  q = model_fn(target, y0 - round(y0));
  
  rr = -0.01:0.025:0.01;
  E = zeros(size(img));
  for r = 1:size(img,1)
    for c = 1:size(img,2)
      candidate = single(extract(img, win, [r c]));
      p = model_fn(candidate, [0 0]);
      E(r,c) = dist_fn(p, q);
    end
  end
end
