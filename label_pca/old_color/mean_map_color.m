function mu_map = mean_map_color(maps, labels)
  %% map to label space
  if ~exist('labels'), labels = single(max(unique(maps{1}))); end

  for i = 1:length(maps)
    lmaps{i} = label_map(maps{i}, labels);
    lmaps{i} = reshape(lmaps{i}, [size(maps{i}) labels-1]);
  end

  %% compute mean map
  mu = m_inv_c(lmaps{1}, @label_unmap, labels);
  for i = 2:length(lmaps)
    mu(:,:,:,2) = m_inv_c(lmaps{i}, @label_unmap, labels);
    mu = sum(mu, 4);
  end
  mu_map = uint8(mu / length(lmaps));

  %% compute mean label map
%   mu = lmaps{1};
%   for i = 2:length(maps)
%     mu(:,:,:,2) = lmaps{i};
%     mu = sum(mu, 4);
%   end
%   mu_lmap = mu / length(maps);


end


