function E = landscape(img, win, y0, model_fn, dist_fn)
  target = single(extract(img, win, round(y0)));
  q = model_fn(target, y0 - round(y0));
  
  E = zeros(size(img));
  for r = 1:size(img,1)
    for c = 1:size(img,2)
      candidate = single(extract(img, win, [r c]));
      p = model_fn(candidate, [0 0]);
      E(r,c) = dist_fn(p, q);
    end
  end
