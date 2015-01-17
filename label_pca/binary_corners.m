function C = binary_corners(label_cnt, cl)
  if ~exist('cl'), cl = 'double'; end
  C = zeros(label_cnt, label_cnt-1, cl);
  C(2:label_cnt+1:end) = 1;
end
