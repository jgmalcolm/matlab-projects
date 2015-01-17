function fn = filter_ma(n)
% FILTER_MA Moving average filter closure (ASSUME: 1D or 2D signal)
%
% FN = FILTER_MA(N) Average of last N time samples.
%
% Example:
%  >> fn = filter_ma(3);
%  >> y = fn(4); % 4
%  >> y = fn(5); % 4.5
%  >> y = fn(6); % 5
%  >> y = fn(1); % 4
%  >> y = fn(7); % 4.66
%  >> fn = filter_ma(2);
%  >> y = fn([1 1]'); % [1 1]'
%  >> y = fn([1 3]'); % [1 2]'
%  >> y = fn([3 2]'); % [2 2.5]'
  
  win = []; % window
  i = 1;
  
  inc = @(x) max(mod(x + 1, n+1),1);

  function y = filter(u)
    if isempty(win)
      [y win] = deal(u);
      return
    end
    i = inc(i); % rotate window
    win(:,:,i) = u;
    y = mean(win,ndims(u)+1);
  end
  
  fn = @filter; % return closure reference
end
