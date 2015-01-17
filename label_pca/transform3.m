function m_ = transform3(m, pose)
  x = pose(1); y = pose(2); z = pose(3);
  u = pose(4); v = pose(5); w = pose(6);
  a = pose(7); b = pose(8); c = pose(9);
  p =pose(10); q =pose(11); r =pose(12);
  
  T = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
      [u 0 0 0; 0 v 0 0; 0 0 w 0; 0 0 0 1] * ... % scale
      [1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1] * ... % about x-axis
      [cos(b) 0 sin(b) 0; 0 1 0 0; -sin(b) 0 cos(b) 0; 0 0 0 1] * ... % about y-axis
      [cos(c) -sin(c) 0 0; sin(c) cos(c) 0 0; 0 0 1 0; 0 0 0 1] * ... % about z-axis
      [1 0 0 0; p 1 0 0; p 0 1 0; 0 0 0 1]' * ... % x-shear
      [1 q 0 0; 0 1 0 0; 0 q 1 0; 0 0 0 1]' * ... % y-shear
      [1 0 r 0; 0 1 r 0; 0 0 1 0; 0 0 0 1]';      % z-shear

  % form coordinates for interpolation
  sz = size(m);
  [xx yy zz] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
  coords = [xx(:)'; yy(:)'; zz(:)'];
  coords(4,:) = 1;
  clear xx yy zz;
  
  % backproject coordinates
  for i = 1:3, coords(i,:) = coords(i,:) - sz(i)/2; end % center about origin
  coords = inv(T) * coords;
  for i = 1:3, coords(i,:) = coords(i,:) + sz(i)/2; end % un-center
  
  for i = 1:size(m,4)
    m_(:,i) = interp3(m(:,:,:,i), ...
                      coords(2,:)', coords(1,:)', coords(3,:)', ... % HACK!
                      'linear', 0);
  end
  m_ = reshape(m_, size(m));
end
