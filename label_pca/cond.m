function P = cond(P, P_)
  sz = size(P);
  label_cnt = sz(end);
  
  P = P .* P_; % condition
  
  % re-normalize
  P = reshape(P, [numel(P)/label_cnt label_cnt]);
  w = sum(P,2) + eps;
  for i = 1:label_cnt
    P(:,i) = P(:,i) ./ w;
  end

  P = reshape(P, sz);
end
