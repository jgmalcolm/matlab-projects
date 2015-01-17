function run_register(pose, register_fn)
  
  if ~exist('register_fn');
    register_fn = @register;
  end
  
  % use first map
  maps = load('images/FOUR', 'maps');
  maps = maps.maps;
  m = maps{1};
  clear maps;
  maps{1} = m;
  
  % transform to get second map
  label_cnt = numel(unique(maps{1}));
  lmap = label_map(maps{1}, label_cnt);
  lmap = reshape(lmap', [size(maps{1}) label_cnt-1]);
  maps{2} = f_inv(transform(lmap, pose));
  
  % register
  [maps_ pose] = register_fn(maps);

end
  
function map = f_inv(lmap)
  sz = size(lmap);
  lmap = reshape(lmap, [prod(sz([1:end-1])) sz(end)]);
  lmap = lmap';
  label_cnt = size(lmap,1) + 1;
  map = label_unmap(lmap, label_cnt);
  map = uint8(reshape(map, sz(1:2)));
end
