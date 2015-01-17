function m_ = transform3_rigid(m, pose)
  [x y z s t p] = deal(pose(1), pose(2), pose(3), ...  % x/y/z translation
                       pose(4), ...                    % isotropic scale
                       pose(5), pose(6));              % theta/phi rotation
  
  T = [1 0 0 x; 0 1 0 y; 0 0 1 z; 0 0 0 1] * ... % translation
      [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1] * ... % scale
      [sin(p)        -cos(p)         0      0; ...
       cos(p)*cos(t)  sin(p)*cos(t) -sin(t) 0; ...
       cos(p)*sin(t)  sin(p)*sin(t)  cos(t) 0; ...
       0              0              0      1]; % rotation

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
