function x_ = binary_map(x, label_cnt, cl)
  if ~exist('cl'), cl = 'double'; end
  C = binary_corners(label_cnt);
  x_ = reshape( C(x,:), [size(x) label_cnt-1]);
end
