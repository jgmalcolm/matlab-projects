function run_project
  s = load('~/src/matlab/images/OCTOPUS', 'masks');
  fn = shape_pca(cell2mat(map(@(x) double(x(:)), s.masks)), 'l');

  sz = size(s.masks{1});
  id = 1;

  sp(1,2,1); imagesc(zeros(sz), [0 1]);
  hold on; contour(s.masks{id}, [.5 .5], 'w', 'LineWidth', 2); hold off;
  axis image off;

  m_ = reshape(fn.preimage(double(flat(s.masks{id}))), sz);
  sp(1,2,2); imagesc(m_, [0 1]);
  hold on; contour(m_, [.5 .5], 'w', 'LineWidth', 2); hold off;
  axis image off
