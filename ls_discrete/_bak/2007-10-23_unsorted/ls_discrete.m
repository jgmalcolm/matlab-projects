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


% move curve C according to corresponding speeds S
function [phi C] = evolve(phi, C, S, h)
  %- update phi
  phi(C(S > 0)) = -1;
  phi(C(S < 0)) =  1;

  % shift matrices
  shiftD = @(M) M([1 1:end-1],:);
  shiftL = @(M) M(:,[2:end end]);
  shiftR = @(M) M(:,[1 1:end-1]);
  shiftU = @(M) M([2:end end],:);

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
  
  C = find(phi == 0); % new curve
end
