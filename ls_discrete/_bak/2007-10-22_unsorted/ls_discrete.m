function [phi C] = ls_discrete(phi, C, max_iter, iter_e, iter_s, h)
  if iter_s ~= 0, error('smoothing not yet implemented'); end

  for i = 1:max_iter  % only evolve, no smoothing
    %- dilate
    C = order_curve(size(phi), C);
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


% move curve C according to corresponding speeds S
function [phi C] = evolve(phi, C, S, h)
  %- update phi
  phi(C(S > 0)) = -1;
  phi(C(S < 0)) =  1;

  %- mark as interface off-contour points with neighbors of opposite sign
  Nu = shiftD(phi); Nd = shiftU(phi);
  Nl = shiftR(phi); Nr = shiftL(phi);
  is_opposite = Nu.*phi < 0 | Nd.*phi < 0 | Nl.*phi < 0 | Nr.*phi < 0;
  is_opposite(C) = 0; % ignore contour points
  phi(find(is_opposite)) = 0;  % mark those as interface
  
  %- maintain minimal interface
  Nu = shiftD(phi); Nd = shiftU(phi);
  Nl = shiftR(phi); Nr = shiftL(phi);
  all_pos = phi == 0 & (Nu >= 0 & Nd >= 0 & Nl >= 0 & Nr >= 0);
  all_neg = phi == 0 & (Nu <= 0 & Nd <= 0 & Nl <= 0 & Nr <= 0);
  phi(all_pos) =  1;
  phi(all_neg) = -1;
  
  %- move in and move out
  h.move_in(uint32(C(S > 0)'));
  h.move_out(uint32([C(S < 0)' find(all_pos)']));
  
  C = find(phi == 0);
end

function M = shiftD(M)
  M = M([1 1:end-1],:);
end
function M = shiftL(M)
  M = M(:,[2:end end]);
end
function M = shiftR(M)
  M = M(:,[1 1:end-1]);
end
function M = shiftU(M)
  M = M([2:end end],:);
end


function C = order_curve(sz, C)
  % compute pair-wise L2 squared distance
  [C_(1,:) C_(2,:)] = ind2sub(sz, single(C));
  C2 = sum(C_.^2, 1);
  p = size(C_,2);
  d2 = repmat(C2, p, 1) - 2*C_'*C_+ repmat(C2', 1, p); % expand inner product
  
  C = C(1); % rebuild
  C(1,:) = Inf;
  for i = 1:length(C)-1
    [val ind] = min(d2(:,i));
    C(end+1) = C(ind)
  end
  
end
