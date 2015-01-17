function [phi, C] = mask2phi(mask)
% MASK2PHI Generate phi from mask for ls_switching.
%  mask -- 1 iff inside, 0 iff outside
% Note, interface considered inside.
%
% >> [phi C] = mask2phi(mask);
%
% (See ls_switching)
  
  phi = -(2*mask - 1);  % default: far
  
  % determine neighbors
  Nu = shiftD(phi); Nd = shiftU(phi);
  Nl = shiftR(phi); Nr = shiftL(phi);

  % determine outside points with neighbors inside
  is_near  = phi == -1 & (Nu == 1  | Nd == 1  | Nl == 1  | Nr == 1 );
  
  % mark as zero level set
  phi(is_near) = 0;
  
  % order lists
  C = find(is_near);
  C = get_ordered_list(size(mask), C);
  
  % ensure clockwise parameterization
  [yy xx] = ind2sub(size(phi), C);
  C = single(C);
  area = sum( xx.*yy([2:end 1]) - yy.*xx([2:end 1]) );
  if area > 0
    C = C(end:-1:1); % reparameterize
  end

  % type conversion
  phi = int8(phi);
  C = uint32(C);

end


function [ordered] = get_ordered_list(matsize, unordered)
% Shawn Lankton
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




function shift = shiftL(M)
  shift = [ M(:,2:end) M(:,end) ];
end

function shift = shiftR(M)
  shift = [ M(:,1) M(:,1:end-1) ];
end

function shift = shiftU(M)
  shift = shiftL(M')';
end

function shift = shiftD(M)
  shift = shiftR(M')';
end
