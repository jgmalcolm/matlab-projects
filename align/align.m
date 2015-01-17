function imgs_aligned = align(imgs)
% ALIGN Align a set of 2D binary images via gradient descent on translation,
% rotation, and scale parameters.  Based upon Tsai, et. al. "A Shape-Based
% Approach to the Segmentation of Medical Imagery Using Level Sets".
%
% DATA = ALIGN(DATA) Reads in the images in the specified directory returning
% aligned versions and displaying original and aligned.
%
% Example:
%  >> imgs = align(imgs);
  
  % --== Parameters ==--
  smooth_dt = 0.1;
  smooth_iter = 20;
  % maximum value each parameter may increment by
  step_a = 0.05; % x-translation
  step_b = 0.05; % y-translation
  step_h = 0.001;
  step_t = 0.00001;
  % termination parameters
  term_lookback = 3; % compare against energy calculated how long ago?
  term_epsilon = 0.0001; % consider converged if (dE <= eps)
  term_rescale = 0; % number rescales, min size:=(3/4)^term_rescale
  term_maxtime = 70; % maximum number of iterations
  % plotting parameters
  plot_recent = 70; % plot only last N data points



  % --== Construct Set ==--
  img_cnt = size(imgs,3);
  for i = 1:img_cnt
    img = imgs(:,:,i);
    s(i).original = smooth(img, smooth_iter, smooth_dt, @gauss);
    s(i).saved = img;
    % centroid
    [r c] = find(img);
    r_cent = mean(r);
    c_cent = mean(c);
    % initialize parameters
    s(i).a = round(size(img,2)/2) - c_cent; % x-translation
    s(i).b = round(size(img,1)/2) - r_cent; % y-translation
    s(i).h = 1; % scale
    s(i).t = 0*pi/180; % theta (rotation)
                       % recent history of parameter values for plotting
    hist_s(i).a = []; hist_s(i).b = [];
    hist_s(i).h = []; hist_s(i).t = [];
  end
  


  % --== Display Original Set ==--
  clf;
  colormap gray;
  % display initial set
  proj = zeros(size(s(1).original));
  for i = 1:img_cnt
    subplot(3,img_cnt+1,i);
    [M H R] = transform(s(i).a, s(i).b, s(i).h, s(i).t);
    img = transform_img(s(i).original, M, H, R);
    imagesc(img);
    title(sprintf('Image %d', i));
    proj = proj + img;
  end
  % display projected intensity map
  proj = proj/img_cnt;
  subplot(3,img_cnt+1,img_cnt+1);
  imagesc(proj);
  title('Intensity Projection');




  % --== Gradient Descent ==--

  hist_E = [];
  t = 0;
  while ~terminate(hist_E, t)
    fprintf('---=== t=%d ===---\n', t);
    fprintf('     %8s   %8s   %8s   %8s\n', 'a', 'b', 'h', 'theta');
    % for each image, calculate the parameters
    e = 0; % total energy
    for i = 1:img_cnt
      img_i = s(i).original;

      % calculate derivatives of transformed image
      rescale_fx = (3/4)^term_rescale; % rescale factor
      [M H R] = transform(s(i).a, s(i).b, s(i).h*rescale_fx, s(i).t);
      img_i = transform_img(img_i, M, H, R);
      % grab scaled down part (top left)
      img_i = img_i(1:round(end*rescale_fx), ...
                    2:round(end*rescale_fx));
      s(i).transformed = img_i; % save transformed version for displaying
      img_dx = Dx(img_i);
      img_dy = Dy(img_i);
      rows = size(img_i,1);
      cols = size(img_i,2);
      dI = [ img_dx(:)'; img_dy(:)'; zeros(1,rows*cols) ]; % flatten derivs.
      
      % calculate partial of transform
      Ta = transform_a(s(i).a, s(i).b, s(i).h*rescale_fx, s(i).t);
      Tb = transform_b(s(i).a, s(i).b, s(i).h*rescale_fx, s(i).t);
      Th = transform_h(s(i).a, s(i).b, s(i).h*rescale_fx, s(i).t);
      Tt = transform_t(s(i).a, s(i).b, s(i).h*rescale_fx, s(i).t);
      
      % calculate partial of transformed image (flattened into column)
      [x y] = meshgrid(1:rows, 1:cols);
      coords = [ x(:)'; y(:)'; ones(1,rows*cols) ];
      coords_t = inv(M*H*R)*coords; % un-transform original coordinates
      P = [];
      P(:,1) = sum(dI .* (Ta * coords_t),1)'; % a
      P(:,2) = sum(dI .* (Tb * coords_t),1)'; % b
      P(:,3) = sum(dI .* (Th * coords_t),1)'; % h
      P(:,4) = sum(dI .* (Tt * coords_t),1)'; % theta
      
      % for all the other images ...
      img_i = img_i(:)'; % flatten image into row
      dE = zeros(1,4); % energy gradient for each parameter
      for j = 1:img_cnt
        if i == j, continue; end; % skip self
        [M H R] = transform(s(j).a, s(j).b, s(j).h*rescale_fx, s(j).t);
        img_j = transform_img(s(j).original, M, H, R);
        img_j = img_j(1:round(end*rescale_fx), ...
                      2:round(end*rescale_fx));
        img_j = img_j(:)';
        
        % compute integrals of squared difference and summation
        diff = img_i - img_j;
        sqr_diff = sum(diff.*diff);
        summ = img_i + img_j;
        sqr_summ = sum(summ.*summ);
        
        % compute first part of energy functional
        num = 2*(img_i-img_j)*P; % 2*int[(Ii-Ij)*P]
        den = sqr_summ; % int[(Ii+Ij)^2]
        part1 = num/den;
        
        % compute second part of energy functional
        num1 = 2*sqr_diff; % 2*int[(Ii-Ij)^2]
        num2 = (img_i+img_j)*P; % int[(Ii+Ij)*P]
        den = sqr_summ^2; % int[(Ii+Ij)^2]^2
        part2 = -num1*num2/den;
        
        % combine parts to form gradient
        dE = dE + (part1 + part2);
        
        % calculate total energy
        e = e + sqr_diff/sqr_summ;
      end % for other images
      
      % update pose parameters for img_i via gradient
      if i ~= 1 % keep first image fixed
        s(i).a = s(i).a + step_a/rescale_fx*dE(1);
        s(i).b = s(i).b + step_b/rescale_fx*dE(2);
        s(i).h = s(i).h + step_h*dE(3);
        s(i).t = s(i).t - step_t*dE(4);
      end
      fprintf('%2d:  %8.3f   %8.3f   %8.3f   %8.3f\n', ...
              i, s(i).a, s(i).b, s(i).h, s(i).t);
    end % for each image

    % store histories
    [ hist_E hist_s ] = store_history(hist_E, hist_s, e, s);
    display_current(s, t, hist_s, hist_E);
    fprintf('\n  E(%d) = %f\n', t, e);
    fprintf('\n');
    t = t + 1;
  end % gradient descent
  
  % final transformation
  for i = 1:img_cnt
    [M H R] = transform(s(i).a, s(i).b, s(i).h, s(i).t);
    imgs_aligned(:,:,i) = transform_img(s(i).saved, M, H, R);
  end
  



  function [hist_E, hist_s] = store_history(hist_E, hist_s, e, s)
  % STORE_HISTORY Stores the history of energy and parameters trimming each to
  % term_lookback length
    % determine where to start grabbing from in history so as to truncate
    start = 1;
    if length(hist_E) == plot_recent, start = 2; end

    % grab histories
    hist_E = [ hist_E(start:end) e ];
    for i = 1:img_cnt
      hist_s(i).a = [ hist_s(i).a(start:end) s(i).a ];
      hist_s(i).b = [ hist_s(i).b(start:end) s(i).b ];
      hist_s(i).h = [ hist_s(i).h(start:end) s(i).h ];
      hist_s(i).t = [ hist_s(i).t(start:end) s(i).t ];
    end
  end
  

  function r = terminate(hist_E, t)
  % TERMINATE Returns true if should terminate, that is if energy is
  % relatively constant in recent iterations after reducing step sizes several
  % times.  Parameterized to tolerance for "relatively constant",
  % parameterized for how many iterations to look back at, and how many times
  % to decimate the step sizes by before really terminating (see top of file).

    % enough data yet?
    if length(hist_E) <= term_lookback
      r = false;
      return
    end
    
    % maximum difference looking back a few...
    recent_dE = abs(hist_E((end-term_lookback+1):end) - hist_E(end));
    max_recent_dE = abs(max(recent_dE));
    no_change = (max_recent_dE <= term_epsilon || term_maxtime <= t);
    
    if no_change
      if term_rescale == 0
        % really terminate
        r = true;
      else
        % rescale the images back up
        term_rescale = term_rescale - 1;
        fprintf(' ****** RESCALE ****** (%d left)\n', ...
                term_rescale);
        r = false;
      end
    else
      % energy still changing
      r = false;
    end
  end % terminate



  function [M, H, R] = transform(a, b, h, t)
  % TRANSFORM Standard transform matrices
    M = [ 1 0 a; 0 1 b; 0 0 1 ];
    H = [ h 0 0; 0 h 0; 0 0 1 ];
    R = [ cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1 ];
  end % transform
  function T = transform_a(a, b, h, t)
  % TRANSFORM_A Partial of transform matrix with respect to 'a'.
    T = [ 0 0 1; 0 0 0; 0 0 0 ] * ...
        [ h 0 0; 0 h 0; 0 0 1 ] * ...
        [ cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1 ];
  end % transform
  function T = transform_b(a, b, h, t)
  % TRANSFORM_B Partial of transform matrix with respect to 'b'.
    T = [ 0 0 0; 0 0 1; 0 0 0 ] * ...
        [ h 0 0; 0 h 0; 0 0 1 ] * ...
        [ cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1 ];
  end % transform
  function T = transform_h(a, b, h, t)
  % TRANSFORM_H Partial of transform matrix with respect to 'h'.
    T = [ 1 0 a; 0 1 b; 0 0 1 ] * ...
        [ 1 0 0; 0 1 0; 0 0 0 ] * ...
        [ cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1 ];
  end % transform
  function T = transform_t(a, b, h, t)
  % TRANSFORM_T Partial of transform matrix with respect to 'theta'.
    T = [ 1 0 a; 0 1 b; 0 0 1 ] * ...
        [ h 0 0; 0 h 0; 0 0 1 ] * ...
        [ -sin(t) -cos(t) 0; cos(t) -sin(t) 0; 0 0 0 ];
  end % transform
    

  function trans_img = transform_img(img, M, H, R)
  % TRANSFORM_IMG Affine transformation of an image.  Fill with black
  % intensity pixels (0x00).  Use bilinear interpolation.
  %
  % TRANS_IMG = TRANSFORM_IMG(IMG, M, H, R)
    
    T = maketform('affine', M'*R*H); % note, M transposed, order rearranged
    u = [ 1 size(img,2) ];
    v = [ 1 size(img,1) ];
    trans_img = imtransform(img, T, 'bilinear', ...
                            'UData', u, ...
                            'VData', v, ...
                            'XData', u, ...
                            'YData', v, ...
                            'FillValues', 0); % black
  end % transform_img
    
    
  
  

  function display_current(s, t, hist_s, hist_E)
  % DISPLAY_CURRENT Display the current state of the alignment starting on
  % the second row of the figure.
  %
  % DISPLAY_CURRENT(IS, TIME, HIST_IS, HIST_ENERGY) Display the translated set
  % with time and energy histories.
    
    % display current state of set
    proj = zeros(rows,cols);
    for i = 1:img_cnt
      % plot
      subplot(3,img_cnt+1,img_cnt+1+i);
      imagesc(s(i).transformed);
      proj = proj + s(i).transformed; % accumulate projection
      % plot pose parameters
      subplot(3,img_cnt+1,2*(img_cnt+1)+i);
      cla;
      hold on;
      plot(hist_s(i).a, 'r');
      plot(hist_s(i).b, 'b');
      plot(hist_s(i).h, 'g');
      plot(hist_s(i).t*180/pi, 'k');
      hold off;
    end
    % display current projected intensity map
    subplot(3,img_cnt+1,2*(img_cnt+1));
    proj = proj/img_cnt;
    imagesc(proj);
    
    % display total energy
    subplot(3,img_cnt+1,3*(img_cnt+1));
    plot(hist_E);
    title('Total Energy');
    
    drawnow;
  end % display_current
  
end % align
