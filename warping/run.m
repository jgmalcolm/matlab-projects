function run
  iterations = 100;

  paths;
  clf; colormap gray;
  s = load('OCTOPUS', 'masks');
  
%   mask = get_blob_mask(zeros(size(s.masks{1})));

  mask = false(size(s.masks{1}));
  mask(20:40, 20:40) = true;

  mask = s.masks{33};

  %- initialize
  [phi C] = mask2phi(mask);
  h = warping_speed();
  h.init(s.masks);
  %- evolve
  [phi_ C_] = ls_sparse(phi, C, h, iterations);

  %- display
  clf; imagesc(zeros(size(mask)), [0 1]); axis image off;
  hold on;
  contour(phi,  [0 0], 'r', 'LineWidth', 2); % initial
  contour(phi_, [0 0], 'w', 'LineWidth', 2); % final
  hold off;
end
