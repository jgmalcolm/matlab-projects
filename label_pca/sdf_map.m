function x_ = sdf_map(x, label_cnt)
  x_ = [];
  for i = 2:label_cnt % assumes background is #1
    d = mask2dist(x == i);
    x_ = [x_; d(:)];
  end
  x_ = reshape(x_, [size(x) label_cnt-1]);
end
