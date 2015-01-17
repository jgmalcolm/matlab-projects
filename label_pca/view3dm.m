function view3dm(dm)
  colors = {'', 'blue', 'red', 'cyan'};

  for i = 1:size(dm,4)
    sdf = dm(:,:,:,i);
    patch(isosurface(sdf,0), 'FaceColor', colors{i}, 'EdgeColor', 'none');
  end
  
  daspect([1 1 1]);
  axis tight off;
  camlight
  lighting phong;
  view(-90, 0);

end
