function run_probabilities(m, windows, sigma)
  %p.map = @binary_map; p.unmap = @binary_unmap; p.prob = @binary_probability;
  %p.map = @scalar_map; p.unmap = @scalar_unmap;
  %p.map = @s_map; p.unmap = @s_unmap;
  p.map = @label_map; p.unmap = @label_unmap; p.prob = @label_probability;
  %p.map = @tsai_map; p.unmap = @tsai_unmap;

  label_cnt = double(max(m(:)));
  f     = @(m) p.map(m,   label_cnt);
  f_inv = @(m) p.unmap(m, label_cnt);
  
  img{1} = repmat(reshape([255 255 255],[1 1 3]), size(m));
  img{2} = repmat(reshape([255 0 0],[1 1 3]), size(m));
  img{3} = repmat(reshape([0 255 0],[1 1 3]), size(m));
  img{4} = repmat(reshape([0 0 255],[1 1 3]), size(m));

  clf;
  for i = numel(windows)
    m_ = smooth2(f(m), windows(i), sigma);
    
    p = p.prob(m_, label_cnt);
    for j = 1:label_cnt
%       sp(label_cnt, 1, j);
      imagesc(uint8(img{j}.*p(:,:,[j j j]))); axis image off;
%       print('-dpng',   sprintf('figs/smoothing/p_binary/p%d', j));
%       print('-deps2c', sprintf('figs/smoothing/p_binary/p%d', j));
      print('-dpng',   sprintf('figs/smoothing/p_ls/p%d', j));
      print('-deps2c', sprintf('figs/smoothing/p_ls/p%d', j));
    end
  end
end
