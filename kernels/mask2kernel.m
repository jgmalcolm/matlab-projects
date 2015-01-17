function phi = mask2kernel(m)
  phi = bwdist(~m);
  ind = find(phi);
  phi(ind) = phi(ind) - .5;
  phi = phi / sum(phi(:));
end
