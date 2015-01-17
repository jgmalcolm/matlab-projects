function vol = read_HUVA(fn)
  fid = fopen(fn, 'r');
  vol = fread(fid, inf, '*short');
  fclose(fid);
  n = numel(vol)/(256*256);
  vol = reshape(vol, [256 256 n]);
end
