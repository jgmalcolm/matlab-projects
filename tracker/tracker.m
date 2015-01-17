function s = tracker(s, data, i)

  %- initialize on new image
  %c = s.c_predict(s.centroid);
  %s.speed.set_img(extract(data.imgs{i}, s.win, round(c)));
  s.speed.set_img(extract(data.imgs{i}, s.win, round(s.centroid)));
  s.speed.init(s.phi, s.C);

  %- evolve
  [phi C] = ls_discrete(s.phi, s.C, s.iter_e, s.iter_e, s.iter_s, s.speed);

  s.speed.postprocess(phi);

  %- centroid
  st = s.speed.get_state();
  c_win = st.centroid;
  base_full = round(s.centroid) - s.win - 1;
  s.centroid = s.c_filter(base_full + c_win);
  
  %- display results
  c = s.centroid;
  big_phi = ones(size(data.imgs{1}), 'int8');
  big_phi = embed(big_phi, phi, round(c));
  subplot(2,1,1); imagesc(overlay(data.imgs{i}, big_phi<=0)); axis image off;
  drawnow;

end
