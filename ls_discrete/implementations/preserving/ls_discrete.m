function [phi C] = ls_discrete(phi, C, max_iter, iter_e, iter_s, h)
  if iter_s ~= 0, error('smoothing not yet implemented'); end
  for i = 1:max_iter  % only evolve, no smoothing
    %- dilate
    S = h.init_iteration(phi, uint32(C'));
    S(S < 0) = 0; % only allow dilation
    [phi C] = evolve(phi, C, S, h);
    %- contract
    S = h.init_iteration(phi, uint32(C'));
    S(S > 0) = 0; % only allow contraction
    [phi C] = evolve(phi, C, S, h);
  end
  C = uint32(C);
end


% generate list of (four) neighbors of phi(r,c) with opposite sign
function n4 = gather_n4(phi, r, c)
  n4 = [];
  if 1 < r           && phi(r-1,c)*phi(r,c) < 0, n4 = [n4 [r-1 c]']; end
  if r < size(phi,1) && phi(r+1,c)*phi(r,c) < 0, n4 = [n4 [r+1 c]']; end
  if 1 < c           && phi(r,c-1)*phi(r,c) < 0, n4 = [n4 [r c-1]']; end
  if c < size(phi,2) && phi(r,c+1)*phi(r,c) < 0, n4 = [n4 [r c+1]']; end
end


% order points by who is closest to next point
function pts_ = reorder(pts, next_pt)
  pts_ = [];
  while size(pts,2) > 0
    d2 = (pts(1,:) - next_pt(1)).^2 + (pts(2,:) - next_pt(2)).^2;
    [val ind] = min(d2);
    ind_rest = setxor(ind, 1:size(pts,2));
    
    pts_(:,end+1) = pts(:,ind);
    next_pt = pts(:,ind);
    pts = pts(:,ind_rest);
  end
end

% move curve C according to corresponding speeds S
function [phi C] = evolve(phi, C, S, h)
  %- update phi
  phi(C(S > 0)) = -1;
  phi(C(S < 0)) =  1;
  
  C = single(C); % BUG: ind2sub incorrect if using uint32(C)
  C_ = [];
  
  %- 1. move contour in and out
  for i = 1:length(C)
    [r c] = ind2sub(size(phi), C(i));
    if S(i) == 0, C_(:,end+1) = [r c]'; continue; end % no change

    %- gather off-contour points with neighbors of opposite sign
    n4 = gather_n4(phi, r, c);
    if numel(n4) == 0, continue; end % found any?
    if i == 1, prev_C = C(end); else prev_C = C(i-1); end % wrap around
    [prev_C(1) prev_C(2)] = ind2sub(size(phi), prev_C);
    n4 = reorder(n4, prev_C); % order: closest to previous point on curve
    C_ = [C_ n4]; % append

    %- mark those as interface
    phi(sub2ind(size(phi), n4(1,:), n4(2,:))) = 0;
  end
  
  %- 2. maintain minimal interface
  drop_pts = [];
  for i = 1:length(C_)
    % is surrounded by only interface or positive points?
    if n4_all_nonnegative(phi, C_(:,i))
      drop_pts(end+1) = i;
      phi(C_(1,i), C_(2,i)) = 1; % make positive
    elseif n4_all_nonpositive(phi, C_(:,i))
      drop_pts(end+1) = i;
      phi(C_(1,i), C_(2,i)) = -1; % make negative
    end
  end
  if numel(drop_pts)
    drop_C = sub2ind(size(phi), C_(1,drop_pts), C_(2,drop_pts));
    C_ = C_(:,setxor(1:length(C_), drop_pts)); % drop those points
  else
    drop_C = [];
  end

  %- 3. move in/out
  h.move_in(uint32(C(S > 0)'));
  h.move_out(uint32([C(S < 0)' drop_C]));
  
  C = sub2ind(size(phi), C_(1,:), C_(2,:))';
end


function ret = n4_all_nonnegative(phi, C)
  r = C(1); c = C(2);
  ret = false; % default: false
  if 1 < r           && phi(r-1,c) < 0, return; end
  if r < size(phi,1) && phi(r+1,c) < 0, return; end
  if 1 < c           && phi(r,c-1) < 0, return; end
  if c < size(phi,2) && phi(r,c+1) < 0, return; end
  ret = true; % passed check, all either zero or positive
end
function ret = n4_all_nonpositive(phi, C)
  r = C(1); c = C(2);
  ret = false; % default: false
  if 1 < r           && phi(r-1,c) > 0, return; end
  if r < size(phi,1) && phi(r+1,c) > 0, return; end
  if 1 < c           && phi(r,c-1) > 0, return; end
  if c < size(phi,2) && phi(r,c+1) > 0, return; end
  ret = true; % passed check, all either zero or negative
end
