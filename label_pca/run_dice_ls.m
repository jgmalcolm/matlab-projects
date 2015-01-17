function D = run_dice_ls(lmaps)
  
  n = numel(lmaps);
  label_cnt = size(lmaps{1}, ndims(lmaps{1})) + 1;

  for i = 1:n
    ii = setxor(1:n, i);
    % separate
    lm = lmaps{i};
    lmaps_ = {lmaps{ii}};
    % train, project
    fn = shape_pca(cell2mat(map(@flat, lmaps_)));
    lm_ = fn.preimage(flat(lm));
    lm_ = reshape(lm_, size(lm));
    % compare
    D(i,:) = dice_ls(lm, lm_, label_cnt);
  end
end
