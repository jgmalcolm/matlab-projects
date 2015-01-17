function show_slice(img, lm)
  sl = 26;
  label_cnt = size(lm, 4) + 1;
  
  T = @(m) flipdim(permute(m, [3 1 2 4]), 1);

  % transformed
  img  = T(img);
  lm  = T(lm);
  
  % sliced
  img  = img(:,:,26);
  lm = squeeze(lm(:,:,26,:));
  
  o = norm(mean(label_corners(label_cnt)));

  colors = ' rbcg';

  imagesc(img); axis image off;
  hold on;
  for i = 2:label_cnt
    contour(ls2dist(lm, label_cnt, i), [o o], colors(i), 'LineWidth', 2);
  end
  hold off;
end
