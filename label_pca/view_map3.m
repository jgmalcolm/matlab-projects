function view_map3(m, label_cnt)
  colors = {'', 'blue', 'red', 'cyan'};

  for i = 2:label_cnt
    patch(isosurface(m == i,0), 'FaceColor', colors{i}, 'EdgeColor', 'none');
  end
  daspect([1 1 1]);
  axis tight off;
  camlight
  lighting phong;
  view(-90, 0);
end
