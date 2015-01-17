function phi = mask2dist(mask)
  phi = bwdist(mask) - bwdist(~mask) - .5*(~mask) + .5*mask;
end



