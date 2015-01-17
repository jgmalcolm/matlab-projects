function imgs = normalize(imgs, maps, ref_img, ref_map)
% NORMALIZE Subtract off difference to reference mean intensities for each region
  label_cnt = numel(unique(maps{1}));
  
  % prepare reference means
  for i = 1:label_cnt
    mu(i) = mean(flatten(ref_img(ref_map==i)));
  end

  % subtract off differences
  for j = 1:length(maps)
    img = imgs{j};
    map = maps{j};
    for i = 1:label_cnt
      ind = map==i;
      mu_ = mean(flatten(img(ind)));
      img(ind) = img(ind) - (mu_ - mu(i)); % subtract
    end
    imgs{j} = img;
  end
end
