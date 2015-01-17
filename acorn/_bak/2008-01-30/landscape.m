function L = landscape(img, win, cent, model_fn, dist_fn, kern_fn)
  target = single(extract(img, win, cent));
  q = model_fn(target, [0 0], kern_fn);

  L = zeros(size(img));
  for r = 1:size(img,1)
    for c = 1:size(img,2)
      candidate = single(extract(img, win, [r c]));
      p = model_fn(candidate, [0 0], kern_fn);
      L(r,c) = dist_fn(p, q);
    end
  end
