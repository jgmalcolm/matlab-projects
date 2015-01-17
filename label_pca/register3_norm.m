function pose = register3_norm(maps)
% REGISTER Register a set of binary maps via gradient descent
% POSE = REGISTER3(MAPS)
% Energy:   sum_{i=1:N} |M_i - mu|^2 / |M_i + mu|^2
% transform model:  x/y/z translation, isotropic scale, x/y/z-rotation
  paths;
  %%-- parameters
  
  % timesteps
  p.step = [50 .015 .01];
  
  p.max_iterations = 300;
  
  %%-- initialize
  N = numel(maps);
  which_step = 1;
      
  for i = 1:N
    pose{i} = [0 0 0 1 0 0]; % x/y/z-translation, isotropic scale, theta/phi-rotation
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
  coords = grid_coords(single(size(mu)));
   
  %%-- gradient descent
  for t = 1:p.max_iterations
    E = 0; % total energy
 
    %%-- (2) compute pose parameter gradient for each image
    for i = 1:N
      load(sprintf([tempdir 'transformed_%d'],i));
      
      %%-- calculate spatial derivatives of transformed lmap
      m_x = flat(Dy(m_i))';
      m_y = flat(Dx(m_i))';
      m_z = flat(Dz(m_i))';
      dm = [m_x; m_y; m_z];
      dm(4,:) = 0;
      clear m_x m_y m_z
      
      %% calculate pose partials of transformed lmap
      dT = pose_partials(pose{i});

      %% calculate gradient of I_i wrt P_i
      for j = 1:length(dT)
        P(j,:) = sum(dm .* (dT{j} * coords));
      end
      clear dm dT
      
      %% compute contribution of image to gradient
      % int(|Ii - mu|^2)              int(|Ii + mu|^2)
      m_i_minus = flat(m_i - mu);     m_i_plus =  flat(m_i + mu);
      f = sum(m_i_minus.^2);          g = sum(m_i_plus.^2);
      df = P * m_i_minus;             dg = P * m_i_plus;
      % combine
      E_i = f/g;
      dE = (df*g - f*dg)/g^2;
      clear m_i m_i_minus m_i_plus P f df g dg;
      
      %% accumulate energy contributions
      E = E + E_i;
      
      %%-- (3) update pose parameters
      pose{i} = pose{i} + p.step([1 1 1 2 3 3]) .* dE';
      
      fprintf(['%2d:  ' ...
               '(%6.2f,%6.2f,%6.2f)  (%4.2f)  (%5.2f,%5.2f)\n'], ...
              i, ...
              pose{i}(1), pose{i}(2), pose{i}(3), ...
              pose{i}(4), ...
              pose{i}(5), pose{i}(6));
    end
    fprintf('  E=%f  (T=%d)\n', E, t);
    
    %%-- adaptive step size
    if E >= prev_E
      which_step = which_step + 1;
      if which_step == 5
        break % DONE
      end
      p.step = p.step / 2;
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
  x = pose(1); y = pose(2);  z = pose(3);  % x/y/z translation
  s = pose(4);                             % isotropic scale
  t = pose(5); p = pose(6);                % theta/phi rotation
  
  % x-translation
  dT{1} = [0 0 0 1; 0 0 0 0; 0 0 0 0; 0 0 0 0] * ... % translation
          [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1] * ... % scale
          [sin(p)        -cos(p)         0      0; ...
           cos(p)*cos(t)  sin(p)*cos(t) -sin(t) 0; ...
           cos(p)*sin(t)  sin(p)*sin(t)  cos(t) 0; ...
           0              0              0      1]; % rotation
  % y-translation
  dT{2} = [0 0 0 0; 0 0 0 1; 0 0 0 0; 0 0 0 0] * ... % translation
          [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1] * ... % scale
          [sin(p)        -cos(p)         0      0; ...
           cos(p)*cos(t)  sin(p)*cos(t) -sin(t) 0; ...
           cos(p)*sin(t)  sin(p)*sin(t)  cos(t) 0; ...
           0              0              0      1]; % rotation
  % z-translation
  dT{3} = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0] * ... % translation
          [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1] * ... % scale
          [sin(p)        -cos(p)         0      0; ...
           cos(p)*cos(t)  sin(p)*cos(t) -sin(t) 0; ...
           cos(p)*sin(t)  sin(p)*sin(t)  cos(t) 0; ...
           0              0              0      1]; % rotation

      
  % isotropic scale
  dT{4} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0] * ... % scale
          [sin(p)        -cos(p)         0      0; ...
           cos(p)*cos(t)  sin(p)*cos(t) -sin(t) 0; ...
           cos(p)*sin(t)  sin(p)*sin(t)  cos(t) 0; ...
           0              0              0      1]; % rotation

     
  % theta-rotation
  dT{5} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1] * ... % scale
          [ 0              0              0      0; ...
           -cos(p)*sin(t) -sin(p)*cos(t) -cos(t) 0; ...
            cos(p)*cos(t)  sin(p)*cos(t) -sin(t) 0; ...
            0              0              0      1]; % rotation

  % phi-rotation
  dT{6} = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
          [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1] * ... % scale
          [ cos(p)        sin(p)        0 0; ...
           -sin(p)*cos(t) cos(p)*cos(t) 0 0; ...
           -sin(p)*sin(t) cos(p)*sin(t) 0 0; ...
            0             0             0 1]; % rotation
end



function coords = grid_coords(sz)
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


function m = flat(m)
  m = m(:);
end
