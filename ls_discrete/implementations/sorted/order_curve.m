function C_ = order_curve(sz, C)
  %- compute pairwise distances (only above diagonal since symmetric)
  [C_(1,:) C_(2,:)] = ind2sub(sz, C);
  d2 = Inf(length(C)); % default: infinite distance
  for i = 1:length(C)
    for j = i+1:length(C)
      d2(i,j) = (C_(1,i) - C_(1,j))^2 + (C_(2,i) - C_(2,j))^2;
    end
  end

  %- iteratively accumulate
  ind = 1;
  C_ = zeros(size(C));
  C_(ind) = C(ind);
  keyboard
  for i = 2:length(C)
    [val ind] = min(d2(i-1,:));
    if numel(ind) > 1, error('foo'); end
    C_(i) = C(ind);
    d2(:,ind) = Inf; % kill out that index
  end
end
