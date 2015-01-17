function p = binary_probability(x, label_cnt)
  sz = size(x);
  
  n = numel(x)/(label_cnt-1);
  x = reshape(x, [n label_cnt-1]);
  
  p = [ones(n,1) exp(x)];
  Z = sum(p, 2);
  p = p ./ Z(:,ones(1,label_cnt));

  if label_cnt == 2, sz(end+1) = 1; end
  p = reshape(p, [sz(1:end-1) label_cnt]);
end
