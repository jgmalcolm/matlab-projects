function m = m_inv(lmap, unmap, label_cnt)
  %- reshape
  sz = size(lmap);
  if numel(sz)==2, sz(3) = 1; end % scalar is 1 dim which then gets ignored
  lmap = reshape(lmap, [prod(sz(1:end-1)) sz(end)]);

  %- unmap
  m = unmap(lmap, label_cnt);
  
  %- repackage
  m = reshape(m, sz(1:end-1));
end
