function [fn N] = get_HUVA
  basedir = 'images/volumes/cropped';
  d = dir([basedir '/*.map']);

  fn = @get_map;
  N = length(d);

  function m = get_map(id)
    fid = fopen([basedir '/' d(id).name], 'rb');
    m = fread(fid, inf, '*uint8');
    fclose(fid);
    m = reshape(m, [256 256 numel(m)/256/256]);
  end
end
