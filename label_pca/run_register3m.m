function pose_ = run_register3m(pose)
  
  % use first map
  maps = load('images/3D', 'maps');
  maps = maps.maps;
  
  % transform to get second map
  label_cnt = numel(unique(maps{1}));
  lmap = label_map(maps{1}, label_cnt);
  lmap = reshape(lmap', [size(maps{1}) label_cnt-1]);
  maps{2} = f_inv(transform3(lmap, pose));
  
  % register
  pose_ = register3m(@get_map, 2);
  
  function m = get_map(id)
    m = maps{id};
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
