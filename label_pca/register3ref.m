function pose = register3ref(maps, m_ref)
% REGISTER Register a set of maps to a reference map
% POSE = REGISTER3REF(MAPS, REF_MAP)
% Energy:   sum_{i=1:N} |M_i - M_ref|^2
% using translation, scale, rotation, skew)
  paths;
  %%-- parameters
  
  % for simplex representation
  p.multi.step{1} = [0.0010 0.0000003 0.0000002 0.0000001];
  p.multi.step{2} = [0.0008 0.0000002 0.0000002 0.0000001];
  p.multi.step{3} = [0.0004 0.0000001 0.0000001 0.0000001];
  p.multi.step{4} = [0.0001 0.0000001 0.0000001 0.0000001];
  
  p.step = p.multi.step{1};
  p.max_iterations = 70;
  
  %%-- initialize
  label_cnt = size(m_ref, ndims(m_ref)) + 1;
  N = length(maps);
  which_step = 1;
      
  for i = 1:N
    pose{i} = [0 0 0 1 1 1 0 0 0 0 0 0]; % x/y/z-translation, scale, rotation, shear
    m_i = label_map(maps{i}, label_cnt, 'single');
    
    % write out to disk to conserve RAM
    save(sprintf([tempdir 'transformed_%d'],i), 'm_i');
    save(sprintf([tempdir 'original_%d'],i), 'm_i');
    clear m_i;
  end
  clear maps
  
  %% grid coordinates
  coords = grid_coords(pose{1}, single(size(m_ref)), label_cnt);
   
  %%-- gradient descent
  prev_E = NaN;
  for t = 1:p.max_iterations
    E = 0; % total energy
 
    %%-- (2) compute pose parameter gradient for each image
    for i = 1:N
      load(sprintf([tempdir 'transformed_%d'],i));
      
      %%-- calculate spatial derivatives of transformed lmap
      m_x = Dy(m_i); % (function naming convention)
      m_y = Dx(m_i); % (function naming convention)
      m_z = Dz(m_i);
      dm = [m_x(:)'; m_y(:)'; m_z(:)'];
      dm(4,:) = 0;
      clear m_x m_y m_z;
      
      %% calculate pose partials of transformed lmap
      dT = pose_partials(pose{i});

      %% calculate gradient of I_i wrt P_i
      for j = 1:length(dT)
        P(j,:) = sum(dm .* (dT{j} * coords), 1);
      end
      clear dm dT;
      
      %% compute contribution of each image to gradient for each parameter
      m_i_ = m_i(:) - m_ref(:);
      % int(|Ii - m_ref|^2)
      E_i = sum(m_i_.^2);
      % dp int(|Ii - m_ref|^2)
      dE = P * m_i_; % OPT: dropped '2 *' from gradient
      clear m_i m_i_ P;
      
      %% accumulate energy contributions
      E = E + E_i;
      
      %%-- (3) update pose parameters
      pose{i} = pose{i} + p.step([1 1 1 2 2 2 3 3 3 4 4 4]) .* dE';
      
      fprintf(['%2d:  ' ...
               '(%6.2f,%6.2f,%6.2f)  (%4.2f,%4.2f,%4.2f)  ' ...
               '(%5.2f,%5.2f,%5.2f)  (%5.2f,%5.2f,%5.2f)\n'], ...
              i, ...
              pose{i}(1), pose{i}(2), pose{i}(3), ...
              pose{i}(4), pose{i}(5), pose{i}(6), ...
              pose{i}(7), pose{i}(8), pose{i}(9), ...
              pose{i}(10), pose{i}(11), pose{i}(12));
    end
    fprintf('  E=%f  (T=%d)\n', E, t);
    
    %%-- adaptive step size
    if E >= prev_E
      which_step = which_step + 1;
      if which_step == 5
        break;
      end
      p.step = p.multi.step{which_step};
    else
      prev_E = E;
    end

      
    %%-- (4) transform maps
    for i = 1:N
      load(sprintf([tempdir 'original_%d'],i));
      m_i = transform3(m_i, pose{i});
      save(sprintf([tempdir 'transformed_%d'],i), 'm_i');
      clear m_i;
    end
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
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % about z-axis
          [1 0 0 0; p 1 0 0; p 0 1 0; 0 0 0 1]' * ... % x-shear
          [1 q 0 0; 0 1 0 0; 0 q 1 0; 0 0 0 1]' * ... % y-shear
          [1 0 r 0; 0 1 r 0; 0 0 1 0; 0 0 0 1]';      % z-shear
  % y-translation
  dT{2} = [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 0 0 0] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % about z-axis
          [1 0 0 0; p 1 0 0; p 0 1 0; 0 0 0 1]' * ... % x-shear
          [1 q 0 0; 0 1 0 0; 0 q 1 0; 0 0 0 1]' * ... % y-shear
          [1 0 r 0; 0 1 r 0; 0 0 1 0; 0 0 0 1]';      % z-shear
  % z-translation
  dT{3} = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % about z-axis
          [1 0 0 0; p 1 0 0; p 0 1 0; 0 0 0 1]' * ... % x-shear
          [1 q 0 0; 0 1 0 0; 0 q 1 0; 0 0 0 1]' * ... % y-shear
          [1 0 r 0; 0 1 r 0; 0 0 1 0; 0 0 0 1]';      % z-shear

      
  % x-scale
  dT{4} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % about z-axis
          [1 0 0 0; p 1 0 0; p 0 1 0; 0 0 0 1]' * ... % x-shear
          [1 q 0 0; 0 1 0 0; 0 q 1 0; 0 0 0 1]' * ... % y-shear
          [1 0 r 0; 0 1 r 0; 0 0 1 0; 0 0 0 1]';      % z-shear
  % y-scale
  dT{5} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % about z-axis
          [1 0 0 0; p 1 0 0; p 0 1 0; 0 0 0 1]' * ... % x-shear
          [1 q 0 0; 0 1 0 0; 0 q 1 0; 0 0 0 1]' * ... % y-shear
          [1 0 r 0; 0 1 r 0; 0 0 1 0; 0 0 0 1]';      % z-shear
  % z-scale
  dT{6} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % about z-axis
          [1 0 0 0; p 1 0 0; p 0 1 0; 0 0 0 1]' * ... % x-shear
          [1 q 0 0; 0 1 0 0; 0 q 1 0; 0 0 0 1]' * ... % y-shear
          [1 0 r 0; 0 1 r 0; 0 0 1 0; 0 0 0 1]';      % z-shear
     
     
  % x-rotation
  dT{7} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [0 0 0 0; 0 -sin(a) -cos(a) 0; 0 cos(a) -sin(a) 0; 0 0 0 0] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % z-rot
          [1 0 0 0; p 1 0 0; p 0 1 0; 0 0 0 1]' * ... % x-shear
          [1 q 0 0; 0 1 0 0; 0 q 1 0; 0 0 0 1]' * ... % y-shear
          [1 0 r 0; 0 1 r 0; 0 0 1 0; 0 0 0 1]';      % z-shear
  % y-rotation
  dT{8} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [-sin(b) 0 cos(b) 0; 0 0 0 0; -cos(b) 0 -sin(b) 0; 0 0 0 0] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % about z-axis
          [1 0 0 0; p 1 0 0; p 0 1 0; 0 0 0 1]' * ... % x-shear
          [1 q 0 0; 0 1 0 0; 0 q 1 0; 0 0 0 1]' * ... % y-shear
          [1 0 r 0; 0 1 r 0; 0 0 1 0; 0 0 0 1]';      % z-shear
  % z-rotation
  dT{9} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [-sin(c) -cos(c) 0 0; cos(c) -sin(c) 0 0; 0 0 0 0; 0 0 0 0] * ... % about z-axis
          [1 0 0 0; p 1 0 0; p 0 1 0; 0 0 0 1]' * ... % x-shear
          [1 q 0 0; 0 1 0 0; 0 q 1 0; 0 0 0 1]' * ... % y-shear
          [1 0 r 0; 0 1 r 0; 0 0 1 0; 0 0 0 1]';      % z-shear
      
  % x-shear
  dT{10}= [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % about z-axis
          [0 0 0 0; 1 0 0 0; 1 0 0 0; 0 0 0 0]' * ... % x-shear
          [1 q 0 0; 0 1 0 0; 0 q 1 0; 0 0 0 1]' * ... % y-shear
          [1 0 r 0; 0 1 r 0; 0 0 1 0; 0 0 0 1]';      % z-shear
  % y-shear
  dT{11}= [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % about z-axis
          [1 0 0 0; p 1 0 0; p 0 1 0; 0 0 0 1]' * ... % x-shear
          [0 1 0 0; 0 0 0 0; 0 1 0 0; 0 0 0 0]' * ... % y-shear
          [1 0 r 0; 0 1 r 0; 0 0 1 0; 0 0 0 1]';      % z-shear
  % z-shear
  dT{12}= [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % about z-axis
          [1 0 0 0; p 1 0 0; p 0 1 0; 0 0 0 1]' * ... % x-shear
          [1 q 0 0; 0 1 0 0; 0 q 1 0; 0 0 0 1]' * ... % y-shear
          [0 0 1 0; 0 0 1 0; 0 0 0 0; 0 0 0 0]';      % z-shear
end







function coords = grid_coords(pose, sz, label_cnt)
  [xx yy zz] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
  coords = [xx(:)'; yy(:)'; zz(:)'];
  coords(4,:) = 1;
  coords = repmat(coords, [1 label_cnt-1]);
end



function dM = Dx(M)
  shiftL = @(M,i) M(:,[2:end end],:,i);
  shiftR = @(M,i) M(:,[1 1:end-1],:,i);
  
  for i = 1:size(M,4)
    dM(:,:,:,i) = (shiftL(M,i) - shiftR(M,i))/2;
  end
end

function dM = Dy(M)
  shiftD = @(M,i) M([1 1:end-1],:,:,i);
  shiftU = @(M,i) M([2:end end],:,:,i);
  
  for i = 1:size(M,4)
    dM(:,:,:,i) = (shiftU(M,i) - shiftD(M,i))/2;
  end
end

function dM = Dz(M)
  shiftF = @(M,i) M(:,:,[2:end end],i);
  shiftB = @(M,i) M(:,:,[1 1:end-1],i);

  for i = 1:size(M,4)
    dM(:,:,:,i) = (shiftF(M,i) - shiftB(M,i))/2;
  end
end
