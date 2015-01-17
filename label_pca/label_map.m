function x_ = label_map(x, label_cnt, cl)
  if ~exist('cl'), cl = 'double'; end
  C = label_corners(label_cnt, cl);
  x_ = reshape( C(x,:), [size(x) label_cnt-1]);
end
