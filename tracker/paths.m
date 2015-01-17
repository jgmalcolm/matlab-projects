function paths
  
  m_base = '~malcolm/src/matlab/';
  c_base = '~malcolm/src/c/';

  if ~exist('extract'), addpath([m_base 'lib']); end
  if ~exist('Dx'), addpath([m_base 'lib/diffs']); end
  if ~exist('ls_discrete'), addpath([m_base 'ls_discrete/unsorted']); end
  if ~exist('threshold_speed'), addpath([c_base 'ls_discrete/speeds']); end
  if ~exist('filter_k'), addpath([m_base 'filters']); end
