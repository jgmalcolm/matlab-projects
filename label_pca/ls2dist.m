function d2 = ls2dist(m, label_cnt, lbl)
% LS2DIST Distance from specified label

  sz = size(m);
  C = label_corners(label_cnt);
  n = numel(m)/(label_cnt-1);
  m = reshape(m, [n label_cnt-1]);
  
  % L2 distance: |x - c|
  d2 = sqrt(sum( (m - C(lbl*ones(1,n),:)).^2, 2));

  if label_cnt == 2, sz(end+1) = 1; end
  d2 = reshape(d2, sz(1:end-1));
end
