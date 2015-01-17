function fn = filter_k(k)
% FILTER_K Constant gain difference filter.
%
% FN = FILTER_K(K) Weight the difference between current and last.
%
% Example:
%  >> fn = filter_k(.6);
%  >> y = fn(4); % 4
%  >> y = fn(7); % 5.8
  
  s = []; % saved state
  
  function y = filter(u)
    if isempty(s)
      s = u;
    end
    
    y = s + k * (u - s);
    s = y; % save
  end
  
  fn = @filter; % return closure reference
end
