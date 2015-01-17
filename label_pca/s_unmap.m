function x = s_unmap(x_, label_cnt)
  % return label for which each point is closest
  C = s_corners(label_cnt);
  x = which_corner(x_, C);
end
