function D = run_dice_thresh(maps)
  
  n = numel(maps);
  label_cnt = size(maps{1}, ndims(maps{1})) + 1;

  for i = 1:n
    ii = setxor(1:n, i);
    % separate
    m = maps{i};
    maps_ = {maps{ii}};
    % train, project
    fn = shape_pca(cell2mat(map(@flat, maps_)));
    m_ = fn.preimage(flat(m));
    m_ = reshape(m_, size(m));
    % compare
    m = sdf_unmap(m, label_cnt);
    m_ = sdf_unmap(m_, label_cnt);
    D(i,:) = dice(m, m_, label_cnt);
  end
end
