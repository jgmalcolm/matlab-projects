function [maps pose mu_prob reg_imgs] = register3(maps,dir_inter,compute_prob, imgs)
% REGISTER Register a set of maps via gradient descent in label-space
% (translation, scale, rotation)
% Energy:   sum_{i=1:N} int{|M_i - mu|^2} dA
  if(~isempty(imgs))
      register_imgs = 1;
  else
      register_imgs = 0;
      reg_imgs = [];
  end
 
   
  %paths;
  addpath('~/yogesh_pi/phd/labelPCA/label_pca_3D/lib/diffs');
  %%-- parameters
  % BOX
  %for binary map representation
  %p.box.step = [0.00004 0.000000006 0.00000002 0.00000001];
  
  %for simplex representation and no normalizing term in the energy
%   p.box.step{1} = [0.00008 0.000000006 0.00000002 0.00000001];
%   p.box.step{2} = [0.00004 0.000000006 0.00000002 0.00000001];
%   p.box.step{3} = [0.00004 0.000000003 0.00000002 0.00000001];
%   p.box.step{4} = [0.00002 0.000000001 0.00000001 0.00000001];
  
%for normalized energy
  p.box.step{1} = [200 0.02 0.1 0.01];
  p.box.step{2} = [100 0.015 0.05 0.005];
  p.box.step{3} = [50 0.01 0.02 0.005];
  p.box.step{4} = [10 0.005 0.01 0.001];
  
  p.step = p.box.step{1};
  p.max_iterations = 150;
  p.visualize = false;
  which_step = 1;
  
  %%-- initialize
  labels = unique([maps{1}]);
  label_cnt = numel(labels);
  N = length(maps);
      
  for i = 1:N
    %maps{i} = maps{i} ;
    
    pose{i} = [0 0 0 1 1 1 0 0 0 0 0 0]; % x/y/z-translation, x/y/z-scale,
                                         % x/y/z-rotation, x/y/z-shear
    lmap = single(label_map(maps{i}, label_cnt));
    lmaps.original{i} = reshape(lmap', [size(maps{i}) label_cnt-1]); % HACK
    %lmaps.original{i} = single(maps{i});
    m_i = lmaps.original{i};
    
    %save the intermediate (overwritten later on)
    name = sprintf('transformed_%d.mat',i);
    name = [dir_inter name];
    save(name,'m_i');
    
    %save the original, saves RAM memory required
    name = sprintf('original_%d.mat',i);
    name = [dir_inter name];
    save(name,'m_i');
    
    %compute mean
    if(i > 1)
      mu(:,:,:,:,2) = m_i;
      mu = sum(mu,5);
    else
      mu = m_i;
    end
    clear m_i;
    
  end
  mu = single(mu / N);

  clear maps;
  clear lmap;
 
  prev_E = 1;
     
  %% grid coordinates
  coords = single(grid_coords(pose{1}, size(mu), label_cnt));
    
   
  %%-- gradient descent
  for t = 1:p.max_iterations
    %if p.visualize; plot_integrals(lmaps); drawnow; end

   
      
    E = 0; % total energy
 
    %%-- (2) compute pose parameter gradient for each image
    for i = 1:N
      
      %load from disk to save memory
      name = sprintf('transformed_%d.mat',i);
      name = [dir_inter name];
      %name = sprintf('../../../labelPCA/labelRegistration/simplex_rep_intermediate/transformed_%d.mat',i);
      load(name);
     
      
      %%-- calculate spatial derivatives of transformed lmap
      m_x = single(Dy(m_i)); % (function naming convention)
      m_y = single(Dx(m_i)); % (function naming convention)
      m_z = single(Dz(m_i));
      dm = [m_x(:)'; m_y(:)'; m_z(:)'];
      dm(4,:) = 0;
      
      clear m_x;
      clear m_y;
      clear m_z;
      
      %% calculate pose partials of transformed lmap
      dT = pose_partials(pose{i});

      %% calculate gradient of I_i wrt P_i
      for j = 1:length(dT)
        P(j,:) = sum(dm .* (dT{j} * coords), 1);
        fprintf('completed vol %d, dT{%d}\n',i,j);
      end

      
      %pfor(1:length(dT),' P(%d,:) = sum(dm .* (dT{%d} * coords), 1);');
      %fprintf('completed vol %d \n',i);
      
      clear dm;
      
      
      %% compute contribution of each image to gradient for each parameter
      % int(|Ii - mu|^2)
      E_i = sum((m_i(:) - mu(:)).^2);
      denom = sum((m_i(:) + mu(:)).^2);
      % dp int(|Ii - mu|^2)
      dE{i} = 2 * (P/denom) * (m_i(:) - mu(:));
      % second term in normalized with area energy
      dE2 = 2 * (E_i/denom^2) * P * (m_i(:) + mu(:));
      
      dE{i} = dE{i} - dE2;
      %% accumulate energy contributions
      E = E + E_i/denom;
      
      clear m_i;
      
      %%-- (3) update pose parameters
      pose{i} = pose{i} + p.step([1 1 1 2 2 2 3 3 3 4 4 4]) .* dE{i}';
      
      fprintf(['%2d:  (%6.2f,%6.2f,%6.2f)  (%5.3f,%5.3f,%5.3f)  (%5.3f,%5.3f,%5.3f)  ' ...
               '(%5.3f,%5.3f,%5.3f)\n'], i, ...
              pose{i}(1), pose{i}(2), pose{i}(3), ...
              pose{i}(4), pose{i}(5), pose{i}(6), ...
              pose{i}(7), pose{i}(8), pose{i}(9), ...
              pose{i}(10), pose{i}(11), pose{i}(12));
    end
    fprintf('  E=%f  (T=%d)\n', E, t);

    if(t > 1)
          if(E > prev_E)
              for i =1:N
                 pose{i} = pose{i} - p.step([1 1 1 2 2 2 3 3 3 4 4 4]) .* dE{i}';
              end
              which_step = which_step + 1;
              if(which_step == 5);
                break;
              end
              p.step = p.box.step{which_step};
              
          else
            prev_E = E;
          end
    else
          prev_E = E;
    end
      
    %%-- (4) transform maps
    for i = 1:N
      name = sprintf('original_%d.mat',i);
      name = [dir_inter name];
      load(name);
      
      m_i = transform3(m_i, pose{i});
      name = sprintf('transformed_%d.mat',i);
      name = [dir_inter name];
      save(name,'m_i');
    
      %should we register the images as well ?
      if(register_imgs)
          reg_imgs{i} = transform3(imgs{i}, pose{i});
      end
      
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
  
  %MPI_Finalize();
  E_mis = 0;
  mu_map = f_inv(mu);
  %%-- finalize
  for i = 1:N
    name = sprintf('transformed_%d.mat',i);
    name = [dir_inter name];
    load (name);  
    maps{i} = f_inv(m_i);
    E_mis = E_mis + length(find((single(mu_map)-single(maps{i}))~=0));
  end
  
  maps{N+1} = mu_map;
  fprintf('E_mislabel = %d\n',E_mis);
  
  sz = size(mu);
  if(compute_prob)
      prob = label_probability(mu,label_cnt);
      for j = 1:label_cnt
        mu_prob{j} = reshape(prob(j,:),[sz(1) sz(2) sz(3)]);
      end
  else
      mu_prob = [];
  end
  
  pose = [];
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



function plot_integrals(lmaps)
  label_cnt = size(lmaps.original{1}, 4) + 1;
  
  clf; colormap gray;
  plot(lmaps.transformed);

  function plot(lmaps)
    m = f_inv(lmaps{1});
    for i = 2:length(lmaps)
      m(:,:,:,2) = f_inv(lmaps{i});
      m = sum(m, 4);
    end
    imagesc3(m);
  end
end



function map = f_inv(lmap)
  sz = size(lmap);
  if numel(sz)==3, sz(4) = 1; end
  lmap = reshape(lmap, [prod(sz([1:3])) sz(4)]);
  lmap = lmap'; % TODO: rearrange these lines
  label_cnt = size(lmap,1) + 1;
  map = label_unmap(lmap, label_cnt);
  map = uint8(reshape(map, sz(1:3)));
end



function coords = grid_coords(pose, sz, label_cnt)
  [xx yy zz] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
  coords = [xx(:)'; yy(:)'; zz(:)'];
  coords(4,:) = 1;
  coords = repmat(coords, [1 label_cnt-1]);
end
