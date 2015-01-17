function d_ = transform_all(d, pose)
  
  for i = 1:length(pose)
    d_{i} = transform(double(d{i}), pose{i});
  end

end
