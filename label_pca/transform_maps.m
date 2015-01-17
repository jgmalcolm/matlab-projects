function maps_ = transform_all(map, pose)
  
  %% map to label space
  label_cnt = numel(unique(map));
  for i = 1:length(pose)
    lmaps{i} = label_map(map, label_cnt);
    lmaps{i} = reshape(lmaps{i}', [size(map) label_cnt-1]); % HACK
  end

  %% transform in label space
  for i = 1:length(pose)
    lmaps{i} = transform(lmaps{i}, pose{i});
  end
  
  %% reproduce label maps
  for i = 1:length(lmaps)
    maps_{i} = f_inv(lmaps{i});
  end


  function map = f_inv(lmap)
    sz = size(lmap);
    lmap = reshape(lmap, [prod(sz([1:end-1])) sz(end)]);
    lmap = lmap';
    label_cnt = size(lmap,1) + 1;
    map = label_unmap(lmap, label_cnt);
    map = uint8(reshape(map, sz(1:2)));
  end
end


