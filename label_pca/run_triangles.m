function run_triangles(maps, test_ind)
  if ~exist('test_ind'), test_ind = 1; end
  
  %% parameters
  p.map = @label_map; p.unmap = @label_unmap;
  p.mode = 2;
  p.scales = -.5:.5:.5;

  %% initialize
  n = length(maps);
  [r c] = size(maps{1});
  label_cnt = numel(unique([maps{:}])) - 1;
  % map data
  for i = 1:n
    m = maps{i};
    m(m == 4) = 1; % can only display triangles
    x_(:,i) = f(m);
  end

  %% compute modes of variation
  fn = shape_pca(x_);
  for i = 1:length(p.scales)
    [m x_hat] = map_mode(p.mode, p.scales(i));
    x_hats(:,:,i) = x_hat;
  end
  %% image modes of variation
  clf;
  rr = 9:20;
  cc = 12:22;
  i = 1;
  for r = rr
    for c = cc
      subplot(length(rr), length(cc), i);
      pixel = (c-1)*size(maps{1},1)+r;
      x_hat = reshape(x_hats(:,pixel,:), [size(x_hats,1) numel(p.scales)]);
      draw_triangle(x_hat, label_cnt); axis image off;
      set(gca, 'XLim', [0 1], 'YLim', [0 1]);
      title(sprintf('(%d,%d)', r, c));
      i = i + 1;
    end
  end

  function x_ = f(m)
    x_ = p.map(m(:), label_cnt);
    x_ = x_(:);
  end
  function m = f_inv(x_)
    m = p.unmap(x_, label_cnt);
    m = reshape(m, [r c]);
  end
  
  function [m x_hat] = map_mode(mode, scale)
    x_hat = fn.preimage_mode(mode, scale);
    m = f_inv(x_hat);
    x_hat = reshape(x_hat, [label_cnt-1 numel(x_hat)/(label_cnt-1)]);
  end
end
