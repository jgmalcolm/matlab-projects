function paths

  m = '~malcolm/src/matlab';

  if ~exist('weighted_centroid'), addpath([m '/lib']); end
  if ~exist('Dx'), addpath([m '/lib/diffs']); end

  if ~exist('mask2kernel'), addpath([m '/shape_kernels']); end
  if ~exist('filter_ukf'), addpath([m '/filters']); end

end
