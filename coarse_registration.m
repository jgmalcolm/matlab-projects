function xx = coarse_registration(xx)
  
  xx = map(@transpose, xx); % for cov()

  n = numel(xx);
  
  for ref = 1:n
    if isempty(xx{ref}), continue, end
    % compute reference centroid, covariance, stdev
    m_ref = mean(xx{ref});
    c_ref =  cov(xx{ref});
    s_ref =  std(xx{ref});
    
    % compare against rest of images
    for i = ref+1:n
      if isempty(xx{i}), continue, end
      m_i = mean(xx{i});
      c_i =  cov(xx{i});

      % find angle between covariance matrices
      [a d b] = svd(c_ref * c_i');
      rot = a * eye(3) * b';
      theta_x = atan(-rot(2,3)/rot(3,3));
      theta_y = atan(rot(1,3) * cos(theta_x)/rot(3,3));
      theta_z = atan(-rot(1,2)/rot(1,1));
      
      %%-- translation
      t = m_i - m_ref;
      
      %%-- rotation
      R = [1 0 0;0 cos(-theta_x) -sin(-theta_x);0 sin(-theta_x) cos(-theta_x)] * ...
          [cos(-theta_y) 0 sin(-theta_y);0 1 0;-sin(-theta_y) 0 cos(-theta_y)] * ...
          [cos(-theta_z) -sin(-theta_z) 0;sin(-theta_z) cos(-theta_z) 0;0 0 1];

      A = [inv(R) -t']; % inverse transform
      tpts = transform(xx{i}, A, m_i);
      
      %%-- scaling
      m_i = mean(tpts);
      s_i =  std(tpts);
      
      s = diag(sqrt(s_ref ./ s_i));
      A = [s [0 0 0]'];
      tpts = transform(tpts, A, m_i); % respect origin
      
%       fprintf('%2d-%2d, (%5.1f %5.1f %5.1f)  (%5.1f %5.1f %5.1f)  (%0.3f %0.3f %0.3f)\n', ...
%               ref, i, t, [theta_x theta_y theta_z]*180/pi,s);

      xx{i} = tpts;
    end
  end

  xx = map(@transpose, xx); % undo
end

function x_ = transform(x, A, orig)
  for i = 1:3, x(:,i) = x(:,i) - orig(i); end
  x(:,4) = 1;
  x_ = x * A';
  for i = 1:3, x_(:,i) = x_(:,i) + orig(i); end
end
