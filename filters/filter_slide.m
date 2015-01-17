function fn = filter_slide(n)
% FILTER_WIN Sliding sindow.
%
% FN = FILTER_SLIDE(N) Returns closure for window of size N
%
% Example:
%  >> fn = filter_slide(3);
%  >> y = fn(4); % []
%  >> y = fn(5); % []
%  >> y = fn(6); % []
%  >> y = fn(1); % 4
%  >> y = fn(7); % 5
%  >> y = fn(2); % 6
  
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
