function psi = img2psi(imgs, dist)
  label_cnt = numel(dist);
  C = label_corners(label_cnt);
  
  for i = 1:length(imgs)
    img = imgs{i};
    sz = size(img);

    % form probabilistic image
    for j = 1:label_cnt
      P_img(:,j) = dist{j}(img(:) + 1);
    end

    % normalize to produce weights
    w = sum(P_img, 2);
    for j = 2:label_cnt  % don't care about bg (zeros)
      W(:,j) = P_img(:,j) ./ w;
    end

    % label space is weighted combination of vertices
    m = W(:,2) * C(2,:);
    for j = 3:label_cnt
      m = m + W(:,j) * C(j,:);
    end

    psi{i} = reshape(m, [sz label_cnt-1]);
  end
end
