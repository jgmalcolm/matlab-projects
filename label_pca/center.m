function maps = center(maps)
% CENTER Center images via centroid in label space
% Compute centroid on each image plane in label space and take average
  
  %%-- map to label space
  label_cnt = numel(unique(maps{1}));
  for i = 1:length(maps)
    lmap = label_map(maps{i}, label_cnt);
    lmaps{i} = reshape(lmap', [size(maps{i}) label_cnt-1]); % HACK
  end

  %%-- translate so centroid at center
  for i = 1:length(maps)
    centroid = [0 0];
    for j = 1:label_cnt-1
      [yy xx] = 
    centroid(:,j) = 
    end
  end
