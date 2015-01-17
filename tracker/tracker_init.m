function s = tracker_init(data)

  %%-- Parameters --%%
  % window padding
  win_pad = 8;
  % iterations of evolution
  s.iter_e = 10;
  % iterations of smoothing
  s.iter_s = 0;
  % centroid filter
  centroid_weight_of_current = 1;
  % speed function
  %speed = @threshold_speed;
  speed = @centroid_speed;

  [yy xx] = find(data.mask);

  % centroid relative to full image
  cy = round(mean(yy));
  cx = round(mean(xx));
  s.centroid = [mean(yy) mean(xx)];
  s.c_filter = filter_k(centroid_weight_of_current);
  %[s.c_estimate s.c_predict] = kalman(A, B, C
  
  % window extent
  dx_right = abs(cx - max(xx));
  dx_left  = abs(cx - min(xx));
  dy_down  = abs(cy - max(yy));
  dy_up    = abs(cy - min(yy));
  s.win(2) = max([ dx_right dx_left ]) + win_pad;
  s.win(1) = max([ dy_up    dy_down ]) + win_pad;
  
  % window versions of phi, C, centroid
  mask_win = extract(data.mask, s.win, [cy cx]);
  [phi C] = mask2phi(mask_win);
  s.phi = phi;
  s.C = C;
  
  % speed functional
  s.speed = speed();
  s.speed.set_img(extract(data.imgs{1}, s.win, [cy cx]));
  s.speed.init(phi, C);

end
