function d = distributions(imgs, maps)
  label_cnt = numel(unique(maps{1}));
  N = length(maps);
  
  for j = 1:label_cnt
    d{j} = zeros(256,1);
  end

  for i = 1:N
    img = imgs{i};
    for j = 1:label_cnt
      h = imhist(img(maps{i} == j));
      d{j} = ((i-1)*d{j} + h/sum(h))/i; % recursively update
    end
  end
end
