function [maps pose] = register_tsai(maps)
% REGISTER Register a set of maps via gradient descent in label-space
% Energy:   sum_{i != j} int{|I_i - I_j|^2} dA / int{|I_i + I_j|^2} dA
  paths;
  %%-- parameters
  % MGH
  p.mgh.step.translate = 0.00025;
  p.mgh.step.scale     = 0.00000001;
  p.mgh.step.rotate    = 0.0000002;
  p.mgh.step.shear     = 0.000000005;
  % TSAI
  p.tsai.step.translate = 0.2;
  p.tsai.step.scale     = 0.00005;
  p.tsai.step.rotate    = 0.0004;
  p.tsai.step.shear     = 0.000004;
  % BAR
  p.bar.step.translate = 0.085;
  p.bar.step.scale     = 0.0001 * 0;
  p.bar.step.rotate    = 0.01 * 0;
  p.bar.step.shear     = 0.00001 * 0;
  
  p.step = p.bar.step;
  p.max_iterations = 100;
  p.visualize = true;

  %%-- initialize
  label_cnt = numel(unique([maps{:}])) - 1;
  for i = 1:length(maps)
    pose{i}.translate = [0 0];
    pose{i}.scale     = [1 1];
    pose{i}.rotate    = 0;
    pose{i}.shear     = 0;
    for j = 1:label_cnt
      lmaps{j}.original{i} = (maps{i} == j+1);   % skip background (#1)
      lmaps{j}.transformed{i} = lmaps{j}.original{i};
    end
    % last channel is actually label space representation for visualization
    lmap = label_map(maps{i}, label_cnt+1);
    lmaps{label_cnt+1}.original{i} = reshape(lmap', [size(maps{i}) label_cnt]); % HACK
    lmaps{label_cnt+1}.transformed{i} = lmaps{label_cnt+1}.original{i};
  end


  %%-- gradient descent
  for t = 1:p.max_iterations
    if p.visualize; plot_integrals(lmaps); drawnow; end
    
    E = 0; % total energy
    %%-- (1) compute pose parameter gradient for each image
    for i = 2:length(maps) % hold first fixed
      %% calculate pose partials of transformed lmap
      [Tx Ty Tu Tv Tr Tk] = pose_partials(pose{i});

      %% grid coordinates
      coords = grid_coords(pose{i}, size(maps{1}), label_cnt);
        
      dE = zeros([6 1]); % energy gradient
      for k = 1:label_cnt % for each label
        img_i = lmaps{k}.transformed{i};
        %%-- calculate coordinate partials of transformed lmap
        img_x = Dx(img_i);
        img_y = Dy(img_i);
        d_img = [img_x(:)'; img_y(:)'];
        d_img(3,:) = 0;
        
        %% calculate gradient of (I_i) wrt P_i
        P(1,:) = sum(d_img .* (Tx * coords), 1); % x translation
        P(2,:) = sum(d_img .* (Ty * coords), 1); % y translation
        P(3,:) = sum(d_img .* (Tu * coords), 1); % x scale
        P(4,:) = sum(d_img .* (Tv * coords), 1); % y scale
        P(5,:) = sum(d_img .* (Tr * coords), 1); % rotation
        P(6,:) = sum(d_img .* (Tk * coords), 1); % shear
        
        %% compute contribution of each image to gradient for each parameter
        for j = setxor(i, 1:length(maps))
          img_j = lmaps{k}.transformed{j};
          %% compute energy gradient for each parameter
          % int(|Ii - Ij|^2)
          f = sum((img_i(:) - img_j(:)).^2);
          % dp int(|Ii - Ij|^2)
          df = 2 * P * (img_i(:) - img_j(:));
          % int(|Ii + Ij|^2)
          g = sum((img_i(:) + img_j(:)).^2);
          % dp int(|Ii + Ij|^2)
          dg = 2 * P * (img_i(:) + img_j(:));
          %% accumulate energy contributions
          E = E + f/g;
          dE = dE + (df.*g - f.*dg)./g.^2;
        end
      end
      
      %%-- (2) update pose parameters
      pose{i}.translate = pose{i}.translate + p.step.translate * dE([1 2])';
      pose{i}.scale     = pose{i}.scale     + p.step.scale     * dE([3 4])';
      pose{i}.rotate    = pose{i}.rotate    - p.step.rotate    * dE(5);
      pose{i}.shear     = pose{i}.shear     + p.step.shear     * dE(6);
      fprintf('%2d:  %8.3f   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f\n', ...
              i, pose{i}.translate(1), pose{i}.translate(2), ...
              pose{i}.scale(1), pose{i}.scale(2), ...
              pose{i}.rotate, pose{i}.shear);
    end
    fprintf('  E=%f  (T=%d)\n', E, t);

    %%-- (3) transform maps
    for i = 1:length(maps)
      for j = 1:label_cnt
        lmaps{j}.transformed{i} = transform(lmaps{j}.original{i}, pose{i});
      end
      if p.visualize
        lmaps{end}.transformed{i} = transform(lmaps{end}.original{i}, pose{i});
      end
    end
  end
  
  %%-- final transformation
  for i = 1:length(maps)
    maps{i} = f_inv(transform(lmaps{end}.original{i}, pose{i}));
  end
end


function T = transformation(pose)
  x = pose.translate(1);
  y = pose.translate(2);
  u = pose.scale(1);
  v = pose.scale(2);
  t = pose.rotate;
  k = pose.shear;

  T = [ 1 0 x; 0 1 y; 0 0 1 ] * ...
      [ u 0 0; 0 v 0; 0 0 1 ] * ...
      [ cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1 ] * ...
      [ 1 k 0; k 1 0; 0 0 1];
end


function [Tx Ty Tu Tv Tr Tk] = pose_partials(pose)
  x = pose.translate(1);
  y = pose.translate(2);
  u = pose.scale(1);
  v = pose.scale(2);
  t = pose.rotate;
  k = pose.shear;

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
  clf; colormap gray;
  
  plot(lmaps{end}.transformed);
  
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
  lmap = reshape(lmap, [prod(sz([1:end-1])) sz(end)]);
  lmap = lmap';
  label_cnt = size(lmap,1) + 1;
  map = label_unmap(lmap, label_cnt);
  map = uint8(reshape(map, sz(1:2)));
end


function coords = grid_coords(pose, sz, label_cnt)
  [xx yy] = meshgrid((1:sz(2)) - sz(2)/2, (1:sz(1)) - sz(1)/2);
  coords = [xx(:)'; yy(:)'];
  coords(3,:) = 1;
  coords = repmat(coords, [1 label_cnt-1]);
end
