function report_dice(d)
  
  label_cnt = size(d{1}.sdf, 2);

  label =  reshape(cell2mat(map(@(x) x.label,  d)), [3 label_cnt numel(d)]);
  sdf =    reshape(cell2mat(map(@(x) x.sdf,    d)), [3 label_cnt numel(d)]);
  binary = reshape(cell2mat(map(@(x) x.binary, d)), [3 label_cnt numel(d)]);

  sdf_mu = mean(sdf, 3)
  sdf_var = var(sdf, 1, 3)

  binary_mu = mean(binary, 3)
  binary_var = var(binary, 1, 3)

  label_mu = mean(label, 3)
  label_var = var(label, 1, 3)
end
