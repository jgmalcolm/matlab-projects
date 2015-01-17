function paths

  if ~exist('Dx'), addpath('~malcolm/src/matlab/lib/diffs'); end
  if ~exist('mask2kernel'), addpath('~malcolm/src/matlab/shape_kernels'); end
  if ~exist('weighted_centroid'), addpath('~malcolm/src/matlab/lib'); end
  if ~exist('filter_ma'), addpath('~malcolm/src/matlab/filters'); end
