function fn = filter_ra(alpha)
% FILTER_RA Recursive average filter closure
%
% FN = FILTER_RA(ALPHA) Weight new values alpha, old with (1-alpha).
%
% Example:
%  >> fn = filter_ra(.6);
%  >> y = fn(4); % 4
%  >> y = fn(7); % 5.8
  
  s = []; % saved state
  
  function y = filter(u)
    if isempty(s)
      s = u;
    end
    
    y = alpha*u + (1-alpha)*s;
    s = y; % save
  end
  
  fn = @filter; % return closure reference
end
