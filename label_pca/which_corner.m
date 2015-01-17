function x = which_corner(x_, C)
  label_cnt = size(C,1);
  n = numel(x_)/size(C,2);
  x_ = reshape(x_, [n size(C,2)]);
  
  %% L2 distance of each point to each corner
  for lbl = 1:label_cnt
    d2(:,lbl) = sum( (x_ - C(lbl*ones([1 n]),:)).^2, 2);
  end

  %% return whatever label it each point was closest to
  [d2 x] = min(d2,[],2);
end
