function show_me(masks, n)
  for i = 1:length(masks)
    m = masks{i} == n(1);
    for j = 2:numel(n)
      m = m | (masks{i} == n(j));
    end
    imagesc(m); axis image off; pause
  end
end
