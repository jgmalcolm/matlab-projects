function C = tsai_corners(label_cnt)
  C = zeros([label_cnt label_cnt-1]);
  for i = 2:label_cnt
    C(i,i-1) = 1;
  end
end
