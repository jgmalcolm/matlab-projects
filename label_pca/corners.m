function C = corners(label_cnt, cl)
  if ~exist('cl'), cl = 'double'; end
  C = zeros([label_cnt label_cnt-1], cl);
  for i = 1:label_cnt-1
    C(i+1,:) = C(i,:);
    if i-1 > 0
      C(i+1,i-1) = C(i,i-1)/i;
      C(i+1,i) = sqrt(1 - (C(i,i-1) - C(i+1,i-1))^2);
    else
      C(i+1,i) = 1;
    end
  end
end
