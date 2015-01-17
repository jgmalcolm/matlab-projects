function lm = label_restrict(lm, label_cnt)

  sz = size(lm);
  n = numel(lm)/(label_cnt-1);
  lm = reshape(lm, [n label_cnt-1]);
  
  % center about origin
  c = mean(label_corners(label_cnt));
  lm = lm - c(ones(1,n),:);
  
  % scale down to inscribed sphere (if outside)
  unit_len = norm(c);
  len = sqrt(sum( lm.^2, 2 )); % L2 norm
  is_out = len > unit_len;  % outside?
  lm(is_out,:) = lm(is_out,:) ./ len(is_out,ones(1,label_cnt-1)) * unit_len;
  
  % uncenter
  lm = lm + c(ones(1,n),:);
  
  lm = reshape(lm, sz);
end
