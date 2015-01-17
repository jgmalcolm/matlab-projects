function d = maps2disk(maps)
  d = tempname;
  mkdir(d)

  top = double(max(max([maps{:}])));
  for i = 1:length(maps)
    m = uint8(double(maps{i}) * 255/top);
    imwrite(m, sprintf([d '/map%02d.png'],i));
  end
end
