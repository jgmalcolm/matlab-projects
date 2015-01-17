function [maps pose] = register_fixed(maps, fixed)
% REGISTER Register a set of maps via gradient descent in label-space
% Energy:   sum_{i=1:N} int{|M_i - M_ref|^2} dA
%  where M_ref is one of the maps held fixed as a reference
  paths;
  %%-- parameters
  % MGH
  p.mgh.step = [0*0.0004 0*0.00000001 0.0000008 0*0.000000005];
  % TSAI
  p.tsai.step = [0.002 0.000001 -0.000002 0.00000002];
  % FOUR
  p.four.step = [0.01 0.00005 0.00001 0.00001];
  % BAR
  p.bar.step = [0.01 0 0 0];
  
  p.step = p.mgh.step;
  p.max_iterations = 100;
  p.visualize = true;

  %%-- initialize
  % move fixed map to front
  m = maps{1};
  maps{1} = maps{fixed};
  maps{fixed} = m;
  label_cnt = single(max(unique([maps{:}])));
  for i = 1:length(maps)
    pose{i} = [0 0 1 1 0 0]; % x/y-translation, x/y-scale, rotation, shear
    lmap = label_map(maps{i}, label_cnt);
    lmaps.original{i} = reshape(lmap, [size(maps{i}) label_cnt-1]);
    lmaps.transformed{i} = lmaps.original{i};
  end

  %%-- gradient descent
  for t = 1:p.max_iterations
    if p.visualize; plot_integrals(lmaps); drawnow; end

    %%-- (1) prepare fixed map
    m_fixed = lmaps.transformed{1};

    E = 0; % total energy
    %%-- (2) compute pose parameter gradient for each image
    for i = 2:length(maps)
      m_i = lmaps.transformed{i};
      %%-- calculate coordinate partials of transformed lmap
      m_x = Dx(m_i);
      m_y = Dy(m_i);
      dm = [m_x(:)'; m_y(:)'];
      dm(3,:) = 0;
      
      %% calculate pose partials of transformed lmap
      [Tx Ty Tu Tv Tr Tk] = pose_partials(pose{i});

      %% grid coordinates
      coords = grid_coords(size(maps{1}), label_cnt);
      
      %% calculate gradient of (I_i) wrt P_i
      P(1,:) = sum(dm .* (Tx * coords), 1); % x translation
      P(2,:) = sum(dm .* (Ty * coords), 1); % y translation
      P(3,:) = sum(dm .* (Tu * coords), 1); % x scale
      P(4,:) = sum(dm .* (Tv * coords), 1); % y scale
      P(5,:) = sum(dm .* (Tr * coords), 1); % rotation
      P(6,:) = sum(dm .* (Tk * coords), 1); % shear
      
      %% compute contribution of each image to gradient for each parameter
      % int(|m_i - m_fixed|^2)
      E_i = sum((m_i(:) - m_fixed(:)).^2);
      % dp int(|m_i - m_fixed|^2)
      dE = 2 * P * (m_i(:) - m_fixed(:));
      %% accumulate energy contributions
      E = E + E_i;
      
      %%-- (3) update pose parameters
      pose{i} = pose{i} + p.step([1 1 2 2 3 4]) .* dE';
      fprintf('%2d:  %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f\n', ...
              i, pose{i}(1), pose{i}(2), pose{i}(3), ...
              pose{i}(4), pose{i}(5), pose{i}(6));
    end
    fprintf('  E=%f  (T=%d)\n', E, t);

    %%-- (4) transform maps
    for i = 2:length(maps)
      lmaps.transformed{i} = transform(lmaps.original{i}, pose{i});
    end
  end
  
  %%-- finalize
  for i = 1:length(maps)
    maps{i} = f_inv(lmaps.transformed{i}, label_cnt);
  end
end



function [Tx Ty Tu Tv Tr Tk] = pose_partials(pose)
  x = pose(1); y = pose(2); u = pose(3);
  v = pose(4); t = pose(5); k = pose(6);

  Tx = [0 0 1; 0 0 0; 0 0 0] * ...
       [u 0 0; 0 v 0; 0 0 1] * ...
       [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1] * ...
       [1 k 0; k 1 0; 0 0 1];
  Ty = [0 0 0; 0 0 1; 0 0 0] * ...
       [u 0 0; 0 v 0; 0 0 1] * ...
       [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1] * ...
       [1 k 0; k 1 0; 0 0 1];
  Tu = [1 0 x; 0 1 y; 0 0 0] * ...
       [1 0 0; 0 0 0; 0 0 0] * ...
       [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1] * ...
       [1 k 0; k 1 0; 0 0 1];
  Tv = [1 0 x; 0 1 y; 0 0 0] * ...
       [0 0 0; 0 1 0; 0 0 0] * ...
       [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1] * ...
       [1 k 0; k 1 0; 0 0 1];
  Tr = [1 0 x; 0 1 y; 0 0 0] * ...
       [u 0 0; 0 v 0; 0 0 1] * ...
       [-sin(t) -cos(t) 0; cos(t) -sin(t) 0; 0 0 0] * ...
       [1 k 0; k 1 0; 0 0 1];
  Tk = [1 0 x; 0 1 y; 0 0 0] * ...
       [u 0 0; 0 v 0; 0 0 1] * ...
       [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1] * ...
       [0 1 0; 1 0 0; 0 0 0];
end



function plot_integrals(lmaps)
  label_cnt = size(lmaps.original{1}, 3) + 1;
  
  clf; colormap gray;
  plot(lmaps.transformed);

  function plot(lmaps)
    m = f_inv(lmaps{1}, label_cnt);
    for i = 2:length(lmaps)
      m = m + f_inv(lmaps{i}, label_cnt);
    end
    imagesc(m); axis image off;
  end
end



function map = f_inv(lmap, label_cnt)
  sz = size(lmap);
  lmap = reshape(lmap, [prod(sz([1:end-1])) sz(end)]);
  map = label_unmap(lmap, label_cnt);
  map = uint8(reshape(map, sz(1:2)));
end



function coords = grid_coords(sz, label_cnt)
  [xx yy] = meshgrid(1:sz(2), 1:sz(1));
  coords = [xx(:)'; yy(:)'];
  coords(3,:) = 1;
  coords = repmat(coords, [1 label_cnt-1]);
end
