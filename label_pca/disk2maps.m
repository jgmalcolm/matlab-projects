function maps = disk2maps(d)
  files = dir([d '/map*.png']);

  for i = 1:numel(files)
    m = imread(sprintf([d '/map%02d.png'], i));
    maps{i} = shift_labels_down(m);
  end
end
