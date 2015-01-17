function [maps pose] = register(maps)
% REGISTER Register a set of maps via gradient descent in label-space
% Energy:   sum_{i=1:N} int{|M_i - mu|^2} dA
  paths;
  %%-- parameters
  % PLANES
  p.planes.step = [2e-2 0*2e-5 0*3e-6 0];
  % MGH
  p.mgh.step = [0*0.0004 0*0.00000001 0.0000004 0*0.000000005];
  % TSAI
  p.tsai.step = [0.002 0.000001 -0.000002 0.00000002];
  % FOUR
  p.four.step = [0.01 0.00005 0.00001 0.00001];
  % BAR
  p.bar.step = [0.01 0 0 0];
  % HUVA
  p.huva.step = [5e-3 1e-7 1e-7 0*1e-8];
  
  p.step = p.huva.step;
  p.max_iterations = 30;
  p.visualize = false;
  p.eps = 1/20; % ... of a percent

  %%-- initialize
  label_cnt = single(max(unique([maps{:}])));
  for i = 1:length(maps)
    pose{i} = [0 0 1 1 0 0]; % x/y-translation, x/y-scale, rotation, shear
    lmap = label_map(maps{i}, label_cnt);
    lmaps.original{i} = reshape(lmap, [size(maps{i}) label_cnt-1]);
    lmaps.transformed{i} = lmaps.original{i};
  end

  %%-- gradient descent
  E_ = NaN;
  for t = 1:p.max_iterations
    if p.visualize; plot_integrals(lmaps); drawnow; end

    %%-- (1) compute intrinsic mean
    mu = lmaps.transformed{1};
    for i = 2:length(maps)
      mu(:,:,:,2) = lmaps.transformed{i};
      mu = sum(mu,4);
    end
    mu = mu / length(maps);

    E = 0; % total energy
    %%-- (2) compute pose parameter gradient for each image
    for i = 1:length(maps)
      m_i = lmaps.transformed{i};
      %%-- calculate coordinate partials of transformed lmap
      m_x = Dx(m_i);
      m_y = Dy(m_i);
      dm = [m_x(:)'; m_y(:)'];
      dm(3,:) = 0;
      
      %% calculate pose partials of transformed lmap
      [Tx Ty Tu Tv Tr Tk] = pose_partials(pose{i});

      %% grid coordinates
      coords = grid_coords(size(maps{1}), label_cnt); % TODO: pull out of loop
      
      %% calculate gradient of (I_i) wrt P_i
      P(1,:) = sum(dm .* (Tx * coords)); % x translation
      P(2,:) = sum(dm .* (Ty * coords)); % y translation
      P(3,:) = sum(dm .* (Tu * coords)); % x scale
      P(4,:) = sum(dm .* (Tv * coords)); % y scale
      P(5,:) = sum(dm .* (Tr * coords)); % rotation
      P(6,:) = sum(dm .* (Tk * coords)); % shear
      
      %% compute contribution of each image to gradient for each parameter
      % int(|Ii - mu|^2)
      E_i = sum((m_i(:) - mu(:)).^2);
      % dp int(|Ii - mu|^2)
      dE = 2 * P * (m_i(:) - mu(:));
      %% accumulate energy contributions
      E = E + E_i;
      
      %%-- (3) update pose parameters
      pose{i} = pose{i} + p.step([1 1 2 2 3 4]) .* dE';
      fprintf('%2d: %7.2f  %7.2f  %6.3f  %6.3f  %6.3f  %6.3f\n', i, pose{i});
    end
    fprintf('  E=%f  (T=%d)  (%.2f%%)\n', E, t, 100*(E-E_)/E);
    is_done = abs( 100*(E-E_)/E ) < p.eps;
    E_ = E;

    %%-- (4) transform maps
    for i = 1:length(maps)
      lmaps.transformed{i} = transform(lmaps.original{i}, pose{i});
    end
    %% termination
    if is_done, break; end
  end
  
  %%-- finalize
  for i = 1:length(maps)
    maps{i} = f_inv(lmaps.transformed{i});
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
  keyboard
end



function plot_integrals(lmaps)
  label_cnt = size(lmaps.original{1}, 3) + 1;
  
  clf; colormap gray;
  plot(lmaps.transformed);

  function plot(lmaps)
    m = f_inv(lmaps{1});
    for i = 2:length(lmaps)
      m = m + f_inv(lmaps{i});
    end
    imagesc(m); axis image off;
  end
end



function map = f_inv(lmap)
  sz = size(lmap);
  if numel(sz)==2, sz(3) = 1; end % scalar is 1 dim which then gets ignored
  lmap = reshape(lmap, [prod(sz([1:end-1])) sz(end)]);
  label_cnt = size(lmap,ndims(lmap)) + 1;
  map = label_unmap(lmap, label_cnt);
  map = uint8(reshape(map, sz(1:2)));
end



function coords = grid_coords(sz, label_cnt)
  [xx yy] = meshgrid(1:sz(2), 1:sz(1));
  coords = [xx(:)'; yy(:)'];
  coords(3,:) = 1;
  coords = repmat(coords, [1 label_cnt-1]);
end
