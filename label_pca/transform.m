function m_ = transform(m, pose)
  x = pose(1);
  y = pose(2);
  u = pose(3);
  v = pose(4);
  t = pose(5);
  k = pose(6);
  
  T = [ 1 0 x; 0 1 y; 0 0 1 ]' * ... % HACK: transposed for imtransform--why?
      [ u 0 0; 0 v 0; 0 0 1 ] * ...
      [ cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1 ] * ...
      [1 k 0; k 1 0; 0 0 1];
  T = maketform('affine', T);
  
  sz = size(m);
  xx = [1 sz(2)] - sz(2)/2; % rotate about origin
  yy = [1 sz(1)] - sz(1)/2;
  m_ = imtransform(m, T, 'bilinear', ...
                   'XData', xx, 'YData', yy, ...
                   'UData', xx, 'VData', yy, ...
                   'FillValues', 0);  % fill with background
end
