function x = sdf_unmap(x_, label_cnt)
  sz = size(x_);
  n = numel(x_)/(label_cnt-1);
  x_ = reshape(x_, [n label_cnt-1]);

  % L2 distance of each (negative value) to zero
  for lbl = 2:label_cnt
    d2(:,lbl) = sum( min(x_(:,lbl-1),0).^2, 2);
  end
  % whichever label was most negative
  [d2 x] = max(d2, [], 2);

  if label_cnt == 2, sz(end+1) = 1; end
  x = reshape(uint8(x), sz(1:end-1));
end
