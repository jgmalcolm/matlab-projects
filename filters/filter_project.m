function fn = filter_project(n)
% FILTER_PROJECT Use previous N displacements to project point.
%
% FN = FILTER_PROJECT(N) Average of last N time samples.  Initialize to
% zero displacement.  (default: N=1)
%
% Example:
%  >> fn = filter_project(3);
% >> fn([1 2]'); % [1 2]'
% >> fn([1 3]'); % [1 4]'
% >> fn([1 4]'); % [1 5]'
% >> fn([1 4]'); % [1 4]'
% >> fn([1 4]'); % [1 4]'
  
  if ~exist('n'), n = 1; end
  win = []; % window
  prev_u = []; % previous input
  i = 1;
  
  inc = @(x) max(mod(x + 1, n+1),1);

  function y = filter(u)
    if isempty(win)
      prev_u = u;
    end
    
    win(:,i) = u-prev_u;
    i = inc(i); % rotate window
    y = u + mean(win,2);
    prev_u = u;
  end
  
  fn = @filter; % return closure reference
end
