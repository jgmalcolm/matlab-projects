function d = run_dice(lmaps, psi_test, imgs_test, lmaps_test)
  
  fn = shape_pca(cell2mat(map(@(x) x(:), lmaps)));
  label_cnt = size(lmaps{1}, ndims(lmaps{1}));
  o = norm(mean(label_corners(label_cnt)))^2;
  
  for i = 1:numel(psi_test)
    im = flipdim(permute(imgs_test{i}, [3 1 2]), 1);
    im_slice = im(:,:,25);
    
    clf;
    lm_ground = flipdim(permute(lmaps_test{i}, [3 1 2 4]), 1);
    lm_slice = squeeze(lm_ground(:,:,25,:));
    sp(2,1,1); imagesc(im_slice); hold on;
    contour(ls2iso(lm_slice, 2), [o o], 'r', 'LineWidth', 2);
    contour(ls2iso(lm_slice, 3), [o o], 'b', 'LineWidth', 2);
    title('original');

    lm_alpha = reshape(fn.preimage(flatten(psi_test{i})), size(psi_test{i}));
    lm_alpha = flipdim(permute(lm_alpha, [3 1 2 4]), 1);
    lm_slice = squeeze(lm_alpha(:,:,25,:));
    sp(2,1,2); imagesc(im_slice); hold on;
    contour(ls2iso(lm_slice, 2), [o o], 'r', 'LineWidth', 2);
    contour(ls2iso(lm_slice, 3), [o o], 'b', 'LineWidth', 2);
    title('preimage');

    d(i,:) = dice(lm_ground, lm_alpha);

    print('-dpng', sprintf('figs/final_kpca_Rmu_T%d', i));
  end
end
