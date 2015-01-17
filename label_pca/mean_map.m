function mu = mean_map(maps)
  mu = maps{1};
  for i = 2:length(maps)
    mu = mu + maps{i};
  end
  mu = mu / length(maps);
end


