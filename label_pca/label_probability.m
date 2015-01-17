function p = label_probability(x_, label_cnt)
% P(x==lbl) = exp(-(x - lbl)^2) / sum( exp(-(x-lbl_i)^2) )

  sz = size(x_);
  C = label_corners(label_cnt);
  n = numel(x_)/(label_cnt-1);
  x_ = reshape(x_, [n label_cnt-1]);

  x_ = label_restrict(x_, label_cnt); % restrict to within inscribed sphere

  % L2 distance of each point to each corner
  for lbl = 1:label_cnt
    d2(:,lbl) = sum( (x_ - C(lbl*ones(1,n),:)).^2, 2);
  end
  
  % normalized exponential probabilities
  exp_d2 = exp(-d2 * label_cnt); % TODO: label_cnt used to normalize?
  normalization = sum(exp_d2,2);
  
  p = exp_d2 ./ normalization(:,ones(1,label_cnt));

  if label_cnt == 2, sz(end+1) = 1; end
  p = reshape(p, [sz(1:end-1) label_cnt]);
end
