function fn = filter_init(n)
% FILTER_INIT Fixed average filter closure: average of first N samples.
%
% FN = FILTER_INIT(N)
% 
% Example:
%  >> fn = filter_init(3);
%  >> y = fn(u);
%  >> y = fn(u);
%  >> y = fn(u);
%  >> fn(u); % same as last y

  win = [];
  fn = @filter;

  function y = filter(u)
    if size(win,2) < n
      win(:,end+1) = u(:);
    end
    y = mean(win, 2);
    y = reshape(y, size(u));    
  end
end
