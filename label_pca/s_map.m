function x_ = s_map(x, label_cnt)
  C = s_corners(label_cnt);
  x_ = C(x(:),:);
end
