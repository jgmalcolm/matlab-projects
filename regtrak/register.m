function P = register(img1, img2)
  paths;

  %%-- parameters
  perc_pixels_used = .20;
  subdivisions = 4;
  pyramid_height = 3;
  N_max = 20;
  dP_eps = 0.1;
  
  %%-- (1) construct Gaussian image pyramid
  img{1}{1} = double(img1);
  img{2}{1} = double(img2);
  for i = 2:pyramid_height
    img{1}{i} = reduce(img{1}{i-1});
    img{2}{i} = reduce(img{2}{i-1});
  end

  %%-- (2) find pixels of largest gradient magnitude at coarsest level (SIGM)
  S2 = select_pixels(img{2}{end}, subdivisions, perc_pixels_used);

%   mask = logical(zeros(size(img{2}{end}), 'uint8'));
%   mask(sub2ind(size(img{2}{end}), S2(1,:), S2(2,:))) = 1;
%   img = uint8(img{2}{end});
%   img = overlay(img, mask, 'r');
%   imagesc(img); axis image off;

  P = [1 0 0 0 1 0]'; % initial
  for i = pyramid_height:-1:1
    %%-- (3) estimate motion parameters P
    [P S2] = estimate_motion(P, S2, img{1}{i}, img{2}{i}, N_max, dP_eps);
    %%-- (4) upscale S2
    if i ~= 1  % TODO: necessary since i==1 is last before exiting? 
      S2 = 2*S2-1;
    end
  end

end


function [P S2] = estimate_motion(P, S2, img1, img2, N_max, dP_eps)
  %%-- (1) calculate H from S2
  H = calc_H(S2, img2);
  
  for i = 1:N_max
    %%-- (2) determine S1 from inverse projection of S2
    S1 = project_back(S2, P);

    %%-- (3) keep only valid S1 points (inside image domain)
    ind = 1 <= S1(1,:) & S1(1,:) <= size(img1,1) & ...
          1 <= S1(2,:) & S1(2,:) <= size(img1,2); % valid
    S1 = S1(:,ind); % keep valid
    S2 = S2(:,ind);
    H = H(ind, :);

    %%-- (4) compute I_t
    i1 = interp2(img1, S1(2,:), S1(1,:), 'linear');
    i2 = img2(sub2ind(size(img2), S2(1,:), S2(2,:)));
    I_t = (i1 - i2)';
    %fprintf('[%d/%d] energy: %12.3f\n', i, N_max, sum(abs(I_t)));

    %%-- (5) solve P=inv(H'*H)*H'*I_t
    dP = inv(H'*H)*H'*I_t;

    %%-- (6) update P
    P = P + dP;

    %%-- (7) done if movement below threshold
    if max(abs(dP)) < dP_eps, break, end
  end
  %fprintf('\n');
end

function S1 = project_back(S2, P)
  a = P(1); b = P(2); c = P(3);
  d = P(4); e = P(5); f = P(6);
  S1(2,:) = ( e*S2(2,:) - b*S2(1,:) + (b*f-c*e))/(a*e-b*d);
  S1(1,:) = (-d*S2(2,:) + a*S2(1,:) + (c*d-a*f))/(a*e-b*d);
end



function H = calc_H(S, img)
  [dx dy] = deal(zeros([1 size(S,2)]));
  for i = 1:size(S,2)
    dx(i) = pDx(img, S(1,i), S(2,i));
    dy(i) = pDy(img, S(1,i), S(2,i));
  end
  x = S(2,:); y = S(1,:);
  H = [dx .* x; dx .* y; dx; ...
       dy .* x; dy .* y; dy ]';
end



function reduced = reduce(img)  %-- halve image
  k = [1 5 8 5 1]/20; % approximate gaussian
  img = conv2(img, k' * k, 'same');
  reduced = img(3:2:end-2, 3:2:end-2);
end


function ind = select_pixels(img, subdivisions, perc_pixels_used)
  [r c] = size(img);
  r_h = floor(r/subdivisions);
  c_h = floor(c/subdivisions);
  
  xx = single(Dx(img));
  yy = single(Dy(img));
  mag = uint8(sqrt(xx.^2 + yy.^2));

  ind = [];
  for r = 0:(subdivisions-1)
    for c = 0:(subdivisions-1)
      win = mag(1+r*r_h:r*r_h+r_h, 1+c*c_h:c*c_h+c_h);
      [rr cc] = top_vals(win, perc_pixels_used);
      rr = rr + r*r_h;
      cc = cc + c*c_h;
      ind = [ind [rr cc]'];
    end
  end

  function [rr cc] = top_vals(val, p) %-- indices of top 10% of values
    h = imhist(val);
    H = cumsum(h / sum(h));
    thresh = find(1-p <= H, 1);
    [rr cc] = find(thresh <= val);
  end
end



%-- full derivatives --%
function dx = Dx(d)
  dx = (shiftL(d) - shiftR(d))/2;
end
function dy = Dy(d)
  dy = (shiftU(d) - shiftD(d))/2;
end


%-- shift operations --%
function M_ = shiftD(M)
  M_ = M([1 1:end-1],:);
end
function M_ = shiftL(M)
  M_ = M(:,[2:end end]);
end
function M_ = shiftR(M)
  M_ = M(:,[1 1:end-1]);
end
function M_ = shiftU(M)
  M_ = M([2:end end],:);
end


%-- pointwise derivatives --%
function d = pDx_b(D, r, c)
  d = D(r,c) - D(r,c-1);
end
function d = pDx_f(D, r, c)
  d = D(r,c+1) - D(r,c);
end
function d = pDx(D, r, c)
  if      c == 1,         d = pDx_f(D, r, c);
  elseif  c == size(D,2), d = pDx_b(D, r, c);
  else                    d = (D(r,c+1) - D(r,c-1))/2;
  end
end
function d = pDy_b(D, r, c)
  d = D(r,c) - D(r-1,c);
end
function d = pDy_f(D, r, c)
  d = D(r+1,c) - D(r,c);
end
function d = pDy(D, r, c)
  if      r == 1,         d = pDy_f(D, r, c);
  elseif  r == size(D,1), d = pDy_b(D, r, c);
  else                    d = (D(r+1,c) - D(r-1,c))/2;
  end
end
