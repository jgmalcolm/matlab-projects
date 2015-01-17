function imgs = disk2imgs(d)
  files = dir([d '/img*.png']);

  for i = 1:numel(files)
    imgs{i} = imread(sprintf([d '/img%02d.png'], i));
  end
end
