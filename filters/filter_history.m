function fn = filter_history(n)
% FILTER_HISTORY Window of last few scalar values
%
% FN = FILTER_HISTORY(N) Last N time samples.
%
% Example:
%  >> fn = filter_history(3);
%  >> y = fn(4); % 4
%  >> y = fn(5); % [4 5]
%  >> y = fn(6); % [4 5 6]
%  >> y = fn(1); % [5 6 1]
%  >> y = fn(7); % [6 1 7]
  
  win = []; % window
  
  function y = filter(u)
    if isempty(win)
      [y win] = deal(u);
      return
    end
    win(end+1) = u;
    
    % trim
    if n < length(win)
      win = win(2:n+1);
    end
    y = win;
  end
  
  fn = @filter; % return closure reference
end
