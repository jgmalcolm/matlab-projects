function d = imgs2disk(imgs)
  d = tempname
  mkdir(d)

  for i = 1:length(imgs)
    m = uint8(imgs{i});
    imwrite(m, sprintf([d '/img%02d.png'],i));
  end
end
