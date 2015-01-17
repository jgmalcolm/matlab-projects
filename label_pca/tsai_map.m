function x_ = tsai_map(x, label_cnt)
  C = tsai_corners(label_cnt);
  x_ = C(x(:),:);
end
