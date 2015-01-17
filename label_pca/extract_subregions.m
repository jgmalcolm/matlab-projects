function [imgs_ maps_] = extract_subregions(imgs, maps)
  % determine maximum enclosing window
  c = charlen3(maps2masks(maps)) + 4;  % extra padding
  c(is_odd(c)) = c(is_odd(c)) + 1; % ensure even
  win = c/2;
  
  % grab same size subregion from each img and map
  for i = 1:length(imgs)
    c = centroid3(maps{i} > 1);
    maps_{i} = extract3(maps{i}, win, round(c));
    imgs_{i} = extract3(imgs{i}, win, round(c));
  end
end
