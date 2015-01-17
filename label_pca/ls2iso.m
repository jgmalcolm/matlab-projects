function sdf = ls2iso(m, label_cnt, lbl)
% LS2ISO Isosurface for specified label

  sz = size(m);
  C = label_corners(label_cnt);
  n = numel(m)/(label_cnt-1);
  m = reshape(m, [n label_cnt-1]);
  
  % squared signed distance: |x - c|^2
  sdf = sum( (m - C(lbl*ones(1,n),:)).^2, 2);

  if label_cnt == 2, sz(end+1) = 1; end
  sdf = reshape(sdf, sz(1:end-1));
end
