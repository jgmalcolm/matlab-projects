function P = img2p(img, dist)
  label_cnt = numel(dist);
  C = label_corners(label_cnt);
  
  sz = size(img);

  % form probabilistic image
  for j = 1:label_cnt
    P(:,j) = dist{j}(img(:) + 1);
  end

  % normalize
  w = sum(P, 2) + eps;
  for j = 1:label_cnt
    P(:,j) = P(:,j) ./ w;
  end

  P = reshape(P, [sz label_cnt]);
end
