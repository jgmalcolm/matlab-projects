function m = m_inv_c(lmap, unmap, label_cnt)
  %- unmap
  sz = size(lmap);
  if numel(sz)==2, sz(3) = 1; end % scalar is 1 dim which then gets ignored
  lmap = reshape(lmap, [prod(sz(1:end-1)) sz(end)]);
  m = unmap(lmap, label_cnt);
  
  %- convert to colors
  m = colorize_map(m);
  
  %- repackage
  m = reshape(m, [sz(1:end-1) 3]);
end
