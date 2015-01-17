function x = label_unmap(x_, label_cnt)
  % return label for which each point is closest

  sz = size(x_);
  C = label_corners(label_cnt);
  n = numel(x_)/(label_cnt-1);
  x_ = reshape(x_, [n label_cnt-1]);
  
  %% L2 distance of each point to each corner
  for lbl = 1:label_cnt
    d2(:,lbl) = sum( (x_ - C(lbl*ones(1,n),:)).^2, 2);
  end

  %% return whatever label each point was closest to
  [d2 x] = min(d2, [], 2);
  
  if label_cnt == 2, sz(end+1) = 1; end
  x = reshape(uint8(x), sz(1:end-1));
end
