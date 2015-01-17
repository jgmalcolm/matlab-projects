function view3lm(lm)
  colors = {'', 'blue', 'red', 'cyan'};

  label_cnt = size(lm, 4) + 1;
  o = norm(mean(label_corners(label_cnt)));
  for i = 2:label_cnt
    d = ls2iso(lm, label_cnt, i);
    d = smooth3(d);
    patch(isosurface(d, o), 'FaceColor', colors{i}, 'EdgeColor', 'none');
  end
  
  daspect([1 1 1]);
  axis tight off;
  camlight
  lighting phong;
  view(-90, 0);
end
