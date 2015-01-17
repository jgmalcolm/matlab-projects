function run_fixed(ind)
  
  s = load('images/2d/2D');
  
  for i = ind
    fprintf('----==== fix %d ====----\n', i);
    [maps_{i} pose_{i}] = register_fixed(s.maps, i);
  end
  
  save('images/2D_fixed', 'maps_', 'pose_');
