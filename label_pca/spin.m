function spin
  i = 1;
  el_ = 0;
  for az = -90:5:270
    view(az,el_); drawnow;
    print('-dpng', sprintf('/tmp/out/img_%08d.png', i));
    i = i + 1;
  end
  
  for el = [0:5:90 85:-5:-90 -85:5:0]
    view(az,el); drawnow;
    print('-dpng', sprintf('/tmp/out/img_%08d.png', i));
    i = i + 1;
  end
end    
