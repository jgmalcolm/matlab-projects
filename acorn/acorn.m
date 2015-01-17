function [y1 t] = acorn(img, win, model_fn, weight_fn, dK_fn, y0)
  N_max = 10; % max iterations
  mu = 0.5; % termination threshold: half pixel

  for t = 1:N_max
    %- (1) candidate patch
    p_img = single(extract(img, win, round(y0)));
    p = model_fn(p_img, y0 - round(y0));

    %- (2) kernel derivative
    dK = dK_fn(y0 - round(y0), []); % BUG: should be [Kx Ky]

    %- (3) weights from similarity metric
    w = weight_fn(p_img, p);

    %- (4) new location via weighted centroid
    [cc rr] = meshgrid((y0(2)-win(2)):(y0(2)+win(2)), ...
                       (y0(1)-win(1)):(y0(1)+win(1)));
    y1 = sum([dK.*w.*rr(:) dK.*w.*cc(:)])/(sum(w) + eps);

    if norm(y1 - y0) < mu, break, end
    y0 = y1;
  end

end
