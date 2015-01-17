function [phi C] = ls_discrete(phi, C, max_iter, iter_e, iter_s, h)
  if iter_s ~= 0, error('smoothing not yet implemented'); end
  for i = 1:max_iter  % only evolve, no smoothing
    %- dilate
    C = order_curve_Shawn(size(phi), C);
    S = h.init_iteration(phi, uint32(C'));
    S(S < 0) = 0; % only allow dilation
    [phi C] = evolve(phi, C, S, h);
    %- contract
    C = order_curve_Shawn(size(phi), C);
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
  is_opposite = find(is_opposite);
  phi(is_opposite) = 0;  % mark those as interface
  
  %- maintain minimal interface
  Nu = shiftD(phi); Nd = shiftU(phi);
  Nl = shiftR(phi); Nr = shiftL(phi);
  all_pos = phi == 0 & (Nu >= 0 & Nd >= 0 & Nl >= 0 & Nr >= 0);
  all_neg = phi == 0 & (Nu <= 0 & Nd <= 0 & Nl <= 0 & Nr <= 0);
  phi(all_pos) =  1;
  phi(all_neg) = -1;
  
  %- move in and move out
  h.move_in(uint32(is_opposite)');
  h.move_out(uint32([C(S < 0)' find(all_pos)']));
  
  C = find(phi == 0); % new curve
end



function C_ = order_curve(sz, C)
  %- compute pairwise distances (only above diagonal since symmetric)
  C = single(C);
  [C_(1,:) C_(2,:)] = ind2sub(sz, C);
  d2 = Inf(length(C)); % default: infinite distance
  for i = 1:length(C)
    for j = i+1:length(C)
      d2(i,j) = sqrt((C_(1,i) - C_(1,j))^2 + (C_(2,i) - C_(2,j))^2);
    end
  end
  
  %- iteratively accumulate
  ind = 1;
  C_ = zeros(size(C), 'uint32');
  for i = 1:length(C)
    [val ind] = min(d2(ind,:));
    if numel(ind) > 1, error('foo'); end
    C_(i) = C(ind);
    d2(:,ind) = Inf; % kill out that index
  end
end

function ordered = order_curve_Shawn(matsize, unordered)
  unordered = single(unordered);
  % written by Shawn Lankton
  numpoints = numel(unordered);
  
  [y x] = ind2sub(matsize,unordered);
  distmap = zeros(length(y));
  
  ordered = zeros(size(unordered));
  
  for(i=1:numpoints)
    y1 = y(i);
    x1 = x(i);
    for(j=1:numpoints)
      if(i == j)
        distmap(i,j)=inf;
      else
        y2 = y(j);
        x2 = x(j);
        
        distmap(i,j)=sqrt((x1-x2)^2+(y1-y2)^2);
        distmap(j,i)=distmap(i,j);
      end
    end
  end     

  idx = 1;
  for(i=1:numpoints)
    ordered(i) = unordered(idx);
    candidates = (distmap(idx,:));
    next_idx = find(candidates == min(candidates),1);
    distmap(:,idx) = inf;
    idx = next_idx;
  end
end
