function mask = gen_mask(sz,input)
  
  % Note bene, boundary = [cols; rows]
  boundary = round(init_boundary(round(input)));
  
  mask = ones(sz, 'uint8'); % default: inside
  % mark boundary
  mask(sub2ind(sz, boundary(2,:), boundary(1,:))) = 2;
  % mark outside starting points
  [mask(1,:) mask(end,:) mask(:,1) mask(:,end)] = deal(0);
  
  % propagate inward till encounter boundary
  is_changed = true;
  while is_changed
    [ mask is_changed ] = propagate(mask);
  end

  mask(sub2ind(sz, boundary(2,:), boundary(1,:))) = 1;

  mask = logical(mask);
end
  
  
  
  
  
  
  
  
function [D is_changed] = propagate(D)
% PROPAGATE Assumes edges already classified as outside so no need to check
% boundary conditions.  Returns flag: 1 iff made changes.
  
  Dl = shiftR(D);
  Dr = shiftL(D);
  Du = shiftD(D);
  Dd = shiftU(D);

  % propagate outside if unset and neighbor is set
  is_outside = find(D == 1 & (Dl == 0 | Dr == 0 | Du == 0 | Dd == 0));
  D(is_outside) = 0;
  is_changed = ~isempty(is_outside);

end








function boundary = init_boundary(input)
% INIT_BOUNDARY Determine the initial zero-level set points and from those
% construct an initial signed distance function.
%
% BOUNDARY = INIT_BOUNDARY(INPUT) Uses INPUT as user selected grid points so
% as to not halt for user input thus automating the process.  Proceeds with
% constructing the initial signed distance function.
  
  boundary = [];
  input(:,end+1) = input(:,1);  % make pseudo-cyclic, i.e. first item
                                % duplicated at end

  for i = 1:(size(input,2)-1) % last refences first that's copied to end
    boundary = [ boundary genpts(input(:,i), input(:,i+1)) ];
  end


  function pts = genpts(a, b)
  % GENPTS generates intervening points between a and b
    
    % We want the dimension with largest coordinate difference to be our reference
    % point to control the granularity at which we create the connecting line.
    diffs = b - a;
    max_diff = max(abs(diffs));

    %-- Calculate slope vector m --%
    m = (b - a)/max_diff;
    
    %-- Generate intervening points --%
    pts = [];
    for t = 0:(max_diff-1)
      p = a + m*t;
      pts = [ pts p ];
    end
    if isempty(pts)
      pts = a;
    end
  end % genpts

end % init_boundary






function shift = shiftL(M)
% shiftL(M) shifts the matrix left duplicating the rightmost column
  shift = [ M(:,2:end) M(:,end) ];
end

function shift = shiftR(M)
% shiftR(M) shifts the matrix right duplicating the leftmost column
  shift = [ M(:,1) M(:,1:end-1) ];
end

function shift = shiftU(M)
% shiftU(M) shifts the matrix up duplicating the bottom-most column
  shift = shiftL(M')';
end

function shift = shiftD(M)
% shiftU(M) shifts the matrix up duplicating the bottom-most column
  shift = shiftR(M')';
end
