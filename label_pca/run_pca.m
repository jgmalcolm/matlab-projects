function d = run_pca(dmaps, bmaps, lmaps)
  label_cnt = size(lmaps{1}, ndims(lmaps{1})) + 1;
  
  for test = 1:numel(dmaps)
    test
    d{test}.sdf = run(dmaps, @sdf_unmap, label_cnt, test);
    d{test}.binary = run(bmaps, @binary_unmap, label_cnt, test);
    d{test}.label = run(lmaps, @label_unmap, label_cnt, test);
  end
end

function d = run(maps, fn_unmap, label_cnt, test)
  lm_test  = maps{test};
  m_test  = fn_unmap(lm_test, label_cnt);

  maps = {maps{setxor(test, 1:numel(maps))}};

  fn_lpca = shape_lpca(cell2mat(map(@(x) x(:), maps)));
  lm_test_lpca = reshape(fn_lpca.preimage(lm_test(:)), size(lm_test));
  m_test_lpca = fn_unmap(lm_test_lpca, label_cnt);
  d(1,:) = dice(m_test, m_test_lpca, label_cnt);

  fn_lpca_ = shape_pca(cell2mat(map(@(x) x(:), maps)), 'linear');
  lm_test_lpca_ = reshape(fn_lpca_.preimage(lm_test(:)), size(lm_test));
  m_test_lpca_ = fn_unmap(lm_test_lpca_, label_cnt);
  d(2,:) = dice(m_test, m_test_lpca_, label_cnt);

  fn_kpca = shape_pca(cell2mat(map(@(x) x(:), maps)), 'kernel');
  lm_test_kpca = reshape(fn_kpca.preimage(lm_test(:)), size(lm_test));
  m_test_kpca = fn_unmap(lm_test_kpca, label_cnt);
  d(3,:) = dice(m_test, m_test_kpca, label_cnt);
end
