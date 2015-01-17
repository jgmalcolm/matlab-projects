function pose = register3_binary(maps)
% REGISTER Register a set of maps via gradient descent in binary label space
% POSE = REGISTER3(MAPS)
% Energy:   sum_{i=1:N} |M_i - mu|^2
% transform model:  x/y/z translation, isotropic scale, x/y/z-rotation
  paths;
  %%-- parameters
  
  % timesteps
  p.steps{1} = [0.0010 0.00000015 0.0000002];
  p.steps{2} = [0.0008 0.00000010 0.0000002];
  p.steps{3} = [0.0004 0.00000005 0.0000001];
  p.steps{4} = [0.0002 0.00000005 0.0000001];
  
  p.step = p.steps{1};
  p.max_iterations = 300;
  
  %%-- initialize
  N = numel(maps);
  which_step = 1;
      
  for i = 1:N
    pose{i} = [0 0 0 1 0 0 0]; % x/y/z-translation, isotropic scale, x/y/z-rotation
    m_i = maps{i};
    
    % write out to disk to conserve RAM
    save(sprintf([tempdir 'transformed_%d'],i), 'm_i');
    save(sprintf([tempdir 'original_%d'],i), 'm_i');
    
    %compute mean
    if i == 1
      mu = m_i;
    else
      mu(:,:,:,:,2) = m_i;
      mu = sum(mu,5);
    end
    clear m_i;
    
  end
  mu = mu / N;
  clear maps;
 
  prev_E = NaN;
     
  %% grid coordinates
  coords = grid_coords(pose{1}, single(size(mu)));
   
  %%-- gradient descent
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
      clear dm;
      
      %% compute contribution of each image to gradient for each parameter
      m_i_ = m_i(:) - mu(:);
      % int(|Ii - mu|^2)
      E_i = sum(m_i_.^2);
      % dp int(|Ii - mu|^2)
      dE = P * m_i_;
      clear m_i m_i_ P;
      
      %% accumulate energy contributions
      E = E + E_i;
      
      %%-- (3) update pose parameters
      pose{i} = pose{i} + p.step([1 1 1 2 3 3 3]) .* dE';
      
      fprintf(['%2d:  ' ...
               '(%6.2f,%6.2f,%6.2f)  (%4.2f)  (%5.2f,%5.2f,%5.2f)\n'], ...
              i, ...
              pose{i}(1), pose{i}(2), pose{i}(3), ...
              pose{i}(4), ...
              pose{i}(5), pose{i}(6), pose{i}(7));
    end
    fprintf('  E=%f  (T=%d)\n', E, t);
    
    %%-- adaptive step size
    if E >= prev_E
      which_step = which_step + 1;
      if which_step == 5
        break;
      end
      p.step = p.steps{which_step};
    else
      prev_E = E;
    end

      
    %%-- (4) transform maps
    for i = 1:N
      load(sprintf([tempdir 'original_%d'],i));
      m_i = transform3_rigid(m_i, pose{i});
      save(sprintf([tempdir 'transformed_%d'],i), 'm_i');
    
      %compute the new mean
      if(i > 1)
        mu(:,:,:,:,2) = m_i;
        mu = sum(mu,5);
      else
        mu = m_i;
      end
    end
    mu = mu / N;
    clear m_i;
  end
end



function dT = pose_partials(pose)
  x = pose(1);  y = pose(2);  z = pose(3);  % translation
  s = pose(4);                              % isotropic scale
  a = pose(5);  b = pose(6);  c = pose(7);  % rotation
  
  % x-translation
  dT{1} = [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0] * ... % translation
          [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1];      % about z-axis
  % y-translation
  dT{2} = [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 0 0 0] * ... % translation
          [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1];      % about z-axis
  % z-translation
  dT{3} = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0] * ... % translation
          [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1];      % about z-axis

      
  % isotropic scale
  dT{4} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1];      % about z-axis

     
  % x-rotation
  dT{5} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1] * ... % scale
          [0 0 0 0; 0 -sin(a) -cos(a) 0; 0 cos(a) -sin(a) 0; 0 0 0 0] * ... % x-rot
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ...   % y-rot
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1];        % z-rot
  % y-rotation
  dT{6} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ...   % x-axis
          [-sin(b) 0 cos(b) 0; 0 0 0 0; -cos(b) 0 -sin(b) 0; 0 0 0 0] * ... % y-axis
          [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1];        % z-axis
  % z-rotation
  dT{7} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1] * ... % scale
          [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ...   % x-axis
          [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ...   % y-axis
          [-sin(c) -cos(c) 0 0; cos(c) -sin(c) 0 0; 0 0 0 0; 0 0 0 0];      % z-axis
end



function coords = grid_coords(pose, sz)
  [xx yy zz] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
  coords = [xx(:)'; yy(:)'; zz(:)'];
  coords(4,:) = 1;
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
