function C = s_corners(label_cnt)
  C = zeros([label_cnt 2]);
  for i = 2:label_cnt
    C(i,1) = cos(pi/2 + 2*pi/(label_cnt-1) * (i-2));
    C(i,2) = sin(pi/2 + 2*pi/(label_cnt-1) * (i-2));
  end
end
