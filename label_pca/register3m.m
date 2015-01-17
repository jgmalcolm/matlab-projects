function pose = register3m(get_map, N)
% REGISTER3M Memory-conscious volume registration
% (translation, scale, rotation)
% Energy:   sum_{i=1:N} int{|M_i - mu|^2} dA
  paths;
  %%-- parameters
  p.box.step = [0.004 0.000001 -0.00001 -0.000001];
  
  p.step = p.box.step;
  p.max_iterations = 50;
  p.visualize = true;

  %%-- initialize
  m = get_map(1);
  label_cnt = numel(unique(m));
  sz = size(m);
  clear m; % MEMORY
  pose = cell([1 N]);
  [pose{:}] = deal([0 0 0 1 1 1 0 0 0 0 0 0]); % x/y/z-translation, x/y/z-scale,
                                               % x/y/z-rotation, x/y/z-shear

  %%-- gradient descent
  for t = 1:p.max_iterations
    %%-- (1) compute intrinsic mean
    mu = transform3(get_lmap(1), pose{1});
    for i = 2:N
      mu(:,:,:,:,2) = transform3(get_lmap(i), pose{i});
      mu = sum(mu,5);
    end
    mu = mu / N;

    if p.visualize; plot_overlay(mu); drawnow; end

    E = 0; % total energy
    %%-- (2) compute pose parameter gradient for each image
    for i = 1:N
      m_i = transform3(get_lmap(i), pose{i});
      %%-- calculate spatial derivatives of transformed lmap
      m_x = Dy(m_i); % (function naming convention)
      m_y = Dx(m_i); % (function naming convention)
      m_z = Dz(m_i);
      dm = [m_x(:)'; m_y(:)'; m_z(:)'];
      dm(4,:) = 0;
      
      %% calculate pose partials of transformed lmap
      dT = pose_partials(pose{i});

      %% grid coordinates
      coords = grid_coords;
      
      %% calculate gradient of I_i wrt P_i
      for j = 1:length(dT)
        P(j,:) = sum(dm .* (dT{j} * coords), 1);
      end

      %% compute contribution of each image to gradient for each parameter
      m = m_i(:) - mu(:);
      % int(|Ii - mu|^2)
      E_i = sum(m.^2);
      % dp int(|Ii - mu|^2)
      dE = 2 * P * m;
      %% accumulate energy contributions
      E = E + E_i;
      
      %%-- (3) update pose parameters
      pose{i} = pose{i} + p.step([1 1 1 2 2 2 3 3 3 4 4 4]) .* dE';
      fprintf(['%2d:  (%6.2f,%6.2f,%6.2f)  (%5.3f,%5.3f,%5.3f)  (%5.3f,%5.3f,%5.3f)  ' ...
               '(%5.3f,%5.3f,%5.3f)\n'], i, pose{i});
    end
    fprintf('  E=%f  (T=%d)\n', E, t);
  end

  function m = get_lmap(id)
    m = get_map(id);
    m = label_map(m, label_cnt);
    m = reshape(m', [sz label_cnt-1]); % TODO: rework to eliminate transpose
  end
  
  function map = f_inv(lmap)
    lmap = reshape(lmap, [prod(sz([1:3])) size(lmap,4)]);
    lmap = lmap'; % TODO: rearrange these lines
    label_cnt = size(lmap,1) + 1;
    map = label_unmap(lmap, label_cnt);
    map = reshape(map, sz(1:3));
  end

  function coords = grid_coords
    [xx yy zz] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
    coords = [xx(:)'; yy(:)'; zz(:)'];
    coords(4,:) = 1;
    coords = repmat(coords, [1 label_cnt-1]);
  end

  function plot_overlay(lmap)
    clf; colormap gray; imagesc3(f_inv(lmap));
  end
end



function dT = pose_partials(pose)
  x = pose(1); y = pose(2); z = pose(3);
  u = pose(4); v = pose(5); w = pose(6);
  a = pose(7); b = pose(8); c = pose(9);
  p =pose(10); q =pose(11); r =pose(12);
  
  % x-translation
  dT{1} = [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) sin(a) 0; 0 -sin(a) cos(a) 0; 0 0 0 1] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
          [1 p p 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * ... % x-shear
          [1 0 0 0; q 1 q 0; 0 0 1 0; 0 0 0 1] * ... % y-shear
          [1 0 0 0; 0 1 0 0; r r 1 0; 0 0 0 1];      % z-shear
  % y-translation
  dT{2} = [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 0 0 0] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) sin(a) 0; 0 -sin(a) cos(a) 0; 0 0 0 1] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
          [1 p p 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * ... % x-shear
          [1 0 0 0; q 1 q 0; 0 0 1 0; 0 0 0 1] * ... % y-shear
          [1 0 0 0; 0 1 0 0; r r 1 0; 0 0 0 1];      % z-shear
  % z-translation
  dT{3} = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) sin(a) 0; 0 -sin(a) cos(a) 0; 0 0 0 1] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
          [1 p p 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * ... % x-shear
          [1 0 0 0; q 1 q 0; 0 0 1 0; 0 0 0 1] * ... % y-shear
          [1 0 0 0; 0 1 0 0; r r 1 0; 0 0 0 1];      % z-shear
  % x-scale
  dT{4} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] * ... % scale
          [1 0 0 0; 0 cos(a) sin(a) 0; 0 -sin(a) cos(a) 0; 0 0 0 1] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
          [1 p p 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * ... % x-shear
          [1 0 0 0; q 1 q 0; 0 0 1 0; 0 0 0 1] * ... % y-shear
          [1 0 0 0; 0 1 0 0; r r 1 0; 0 0 0 1];      % z-shear
  % y-scale
  dT{5} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0] * ... % scale
          [1 0 0 0; 0 cos(a) sin(a) 0; 0 -sin(a) cos(a) 0; 0 0 0 1] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
          [1 p p 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * ... % x-shear
          [1 0 0 0; q 1 q 0; 0 0 1 0; 0 0 0 1] * ... % y-shear
          [1 0 0 0; 0 1 0 0; r r 1 0; 0 0 0 1];      % z-shear
  % z-scale
  dT{6} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0] * ... % scale
          [1 0 0 0; 0 cos(a) sin(a) 0; 0 -sin(a) cos(a) 0; 0 0 0 1] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
          [1 p p 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * ... % x-shear
          [1 0 0 0; q 1 q 0; 0 0 1 0; 0 0 0 1] * ... % y-shear
          [1 0 0 0; 0 1 0 0; r r 1 0; 0 0 0 1];      % z-shear
  % x-rotation
  dT{7} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [0 0 0 0; 0 -sin(a) cos(a) 0; 0 -cos(a) -sin(a) 0; 0 0 0 0] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
          [1 p p 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * ... % x-shear
          [1 0 0 0; q 1 q 0; 0 0 1 0; 0 0 0 1] * ... % y-shear
          [1 0 0 0; 0 1 0 0; r r 1 0; 0 0 0 1];      % z-shear
  % y-rotation
  dT{8} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) sin(a) 0; 0 -sin(a) cos(a) 0; 0 0 0 1] * ... % x-rot
          [-sin(b) 0 cos(b) 0; 0 0 0 0; -cos(b) 0 -sin(b) 0; 0 0 0 0] * ... % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
          [1 p p 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * ... % x-shear
          [1 0 0 0; q 1 q 0; 0 0 1 0; 0 0 0 1] * ... % y-shear
          [1 0 0 0; 0 1 0 0; r r 1 0; 0 0 0 1];      % z-shear
  % z-rotation
  dT{9} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) sin(a) 0; 0 -sin(a) cos(a) 0; 0 0 0 1] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
          [-sin(c) cos(c) 0 0; -cos(c) -sin(c) 0 0; 0 0 0 0; 0 0 0 0]* ...% z-rot
          [1 p p 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * ... % x-shear
          [1 0 0 0; q 1 q 0; 0 0 1 0; 0 0 0 1] * ... % y-shear
          [1 0 0 0; 0 1 0 0; r r 1 0; 0 0 0 1];      % z-shear
  % x-shear
  dT{10}= [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) sin(a) 0; 0 -sin(a) cos(a) 0; 0 0 0 1] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
          [0 1 1 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] * ... % x-shear
          [1 0 0 0; q 1 q 0; 0 0 1 0; 0 0 0 1] * ... % y-shear
          [1 0 0 0; 0 1 0 0; r r 1 0; 0 0 0 1];      % z-shear
  % y-shear
  dT{11}= [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) sin(a) 0; 0 -sin(a) cos(a) 0; 0 0 0 1] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
          [1 p p 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * ... % x-shear
          [0 0 0 0; 1 0 1 0; 0 0 0 0; 0 0 0 0] * ... % y-shear
          [1 0 0 0; 0 1 0 0; r r 1 0; 0 0 0 1];      % z-shear
  % z-shear
  dT{12}= [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) sin(a) 0; 0 -sin(a) cos(a) 0; 0 0 0 1] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
          [1 p p 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * ... % x-shear
          [1 0 0 0; q 1 q 0; 0 0 1 0; 0 0 0 1] * ... % y-shear
          [0 0 0 0; 0 0 0 0; 1 1 0 0; 0 0 0 0];      % z-shear
end


function m_ = transform3(m, pose)
  x = pose(1); y = pose(2); z = pose(3);
  u = pose(4); v = pose(5); w = pose(6);
  a = pose(7); b = pose(8); c = pose(9);
  p =pose(10); q =pose(11); r =pose(12);
  
  T = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
      [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
      [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % x-rot
      [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
      [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
      [1 p p 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] * ... % x-shear
      [1 0 0 0; q 1 q 0; 0 0 1 0; 0 0 0 1] * ... % y-shear
      [1 0 0 0; 0 1 0 0; r r 1 0; 0 0 0 1];      % z-shear
  clear x y z u v w a b c p q r; % MEMORY
  T = single(T); % MEMORY

  % form coordinates for interpolation
  sz = size(m);
  [xx yy zz] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
  coords = single([xx(:)'; yy(:)'; zz(:)']); % MEMORY
  clear xx yy zz; % MEMORY
  coords(4,:) = 1;
  for i = 1:3, coords(i,:) = coords(i,:) - sz(i)/2; end % center about origin
  coords_ = inv(T) * coords;
  clear coords; % MEMORY
  for i = 1:3, coords_(i,:) = coords_(i,:) + sz(i)/2; end % un-center

  for i = 1:size(m,4)
    m_(:,i) = interp3(m(:,:,:,i), coords_(1,:), coords_(2,:), coords_(3,:), ...
                      'linear', 0);
  end
  m_ = reshape(m_, size(m));
end
