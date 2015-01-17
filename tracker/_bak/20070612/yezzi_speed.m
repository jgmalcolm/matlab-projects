function h = yezzi_speed(img)
  %%-- parameters --%%
  nu = 20; % histogram smoothing kernel size
  hist_rate = 0; % rate of refreshing stored histograms [1=reset, 0=keep]
  prior_damp = 4; % edge dampening weight on exponential prior
  g_thresh = 1e-4; % below this, don't move contours
  register_iter = 5;  % use odd counts
  
  % required
  h.init_iteration = @init_iteration;
  h.move_in = @move_in;
  h.move_out = @move_out;
  % extensions
  h.init = @init;
  h.postprocess = @postprocess;
  h.predict_centroid = @predict_centroid;
  h.set_img = @set_img;
  h.get_state = @get_state;
  
  s = [];  % dummy for lexical scoping
  img_bak = []; % previous image for updating distribution
  kern = binomfilt(nu); % for smooth histogram sampling


  %%-- register --%%
  % move img_hat, move phi, and hold img so as to minimize correlation
  % energy.  use final curve from previous frame for comparisons.
  function [phi, C] = register(phi, C)
    [C(1,:) C(2,:)] = ind2sub(size(phi), C);
    %%-- center saved phi
    [yy xx] = find(s.phi_last <= 0);
    saved_c = [mean(yy) mean(xx)];
    phi_last = grab_patch(s.phi_last, round(saved_c), ...
                          (size(phi)-1)/2, int8(1));
    s.img_hat = grab_patch(s.img_hat, round(saved_c), ...
                           (size(phi)-1)/2, 0*s.img_hat(1));
    function do_iter
      %%-- calculate group action 'g' --%%
      ind = find(phi_last <= 0);
      img_dx = Dx(s.img_hat); img_dx = img_dx(ind);
      img_dy = Dy(s.img_hat); img_dy = img_dy(ind);
      img_n = img(ind) / sqrt( dot(img(ind),img(ind)) );
      I_hat = s.img_hat(ind) + eps;
      I_hat_dot = dot(I_hat, I_hat);
      img_I_hat = dot(img_n, I_hat);
      dg_cor(2) = dot(img_dx, img_n)/sqrt(I_hat_dot) ...
                  - img_I_hat*dot(img_dx, I_hat)/I_hat_dot^(3/2);
      dg_cor(1) = dot(img_dy, img_n)/sqrt(I_hat_dot) ...
                  - img_I_hat*dot(img_dy, I_hat)/I_hat_dot^(3/2);
      dg = -dg_cor;
      if abs(dg(1)) < g_thresh && abs(dg(2)) < g_thresh % Hack! jgm
        dt = 0; % don't move
      else
        dt = 1 / (max(abs(dg)) + eps); % at most rounds to one pixel
      end
      g = int32(round(single( dt*dg )));
      s.img_hat = translate(s.img_hat, g, 0);
      phi = translate(phi, g, int8(1)); % pad: outside
      phi_last = translate(phi_last, g, int8(1)); % pad: outside
      C = uint32(int32(C) + g(ones(1,length(C)), :)');
    end
    
    E = @(I, I_hat) dot(I,I_hat) / sqrt(dot(I,I) * dot(I_hat,I_hat));

    for t = 1:register_iter  % use odd counts
      do_iter();
    end
    C = sub2ind(size(phi), C(1,:), C(2,:));
  end

  
  function create_prior(phi)
    mask = zeros(size(phi));
    mask(phi<=0) = 0.01;
    mask = bwdist(mask) - 2;
    %s.prior = exp( -(bwdist(phi <= 0)-2) / prior_damp );
    pindx = find(mask > 0);
    nindx = find(mask <=0 );
    mask(pindx) = exp(-mask(pindx)/prior_damp);
    mask(nindx) = 1;
    s.prior = mask;
  end
  function create_P(phi)
    img_int = uint8(img(:) + 1); % quantize
    hist  = wkeep(conv(imhist(img_int(phi <= 0)), kern), 256);
    s.in.hist = hist_rate*hist + (1-hist_rate)*s.in.hist;
    s.in.hist = s.in.hist / sum(s.in.hist);
    hist = wkeep(conv(imhist(img_int(phi >  0)), kern), 256);
    s.out.hist = hist_rate*hist + (1-hist_rate)*s.out.hist;
    s.out.hist = s.out.hist / sum(s.out.hist);
  end

  function postprocess(phi)
    s.phi_last = phi; % for registration
    create_prior(phi);
    create_P(phi);
  end

  %store the predicted centroid position
  function predict_centroid(cent)
      s.use_cent = 1;
      s.predict_cent = cent;   
  end

  function ret = get_state
    ret = s;
  end
  
  % called at the beginning of each new image
  function [phi, C] = init(phi, C)
    if ~isfield(s, 'img_hat'), s.img_hat = img; end
    if ~isfield(s, 'phi_last'), s.phi_last = phi; end
    s.use_cent = 0;
    
%     [phi C] = register(phi, C); % register via Fitts correlation

    img_int = uint8(img(:)); % quantize
    if ~isfield(s, 'in') || ~isfield(s, 'out') || ~isfield(s, 'prior')
      %-- first time initialization
      s.in.hist  = wkeep(conv(imhist(img_int(phi <= 0)), kern), 256);
      s.in.hist = s.in.hist / sum(s.in.hist);
      s.out.hist = wkeep(conv(imhist(img_int(phi >  0)), kern), 256);
      s.out.hist = s.out.hist / sum(s.out.hist);
      postprocess(phi);
      s.img_hat = img;
    end
    
    %-- probabilistic image
    pin = s.in.hist(img_int+1);
    pout = s.out.hist(img_int+1);
    s.img_P = (128 * pin ./ (pin + pout + eps) .* s.prior(:))';
    
    figure(100);imagesc(reshape(s.img_P,size(phi))); colormap gray;

    %-- inside
    idx = find(phi <= 0);
    s.in.cnt = length(idx);
    s.in.mean = mean(s.img_P(idx));
    % centroid of inside
    [x y] = unpack_point(idx);
    s.centroid = [mean(y) mean(x)];
    fprintf('yezzi/init: %f\n', s.centroid);
    
    %-- outside
    idx = find(phi > 0);
    s.out.cnt = length(idx);
    s.out.mean = mean(s.img_P(idx));
%     fprintf('[%f %f]\n', s.in.mean, s.out.mean);
  end


  

  function S = init_iteration(phi, C)
%     figure(100);
%     imagesc(reshape(s.img_P, size(phi))); colormap gray;

    [xx yy] = ind2sub(size(phi), find(phi <= 0));
    c = [mean(yy) mean(xx)];
    fprintf('yezzi/init_iteration: %f\n', c);
    
    %%-- compute H0 gradient --%%
    u = s.in.mean;
    v = s.out.mean;
    %%-- (1) Yezzi on current image
    f_y =  5*(u-v) * ((s.img_P(C)-v)/s.out.cnt + ...
                      (s.img_P(C)-u)/s.in.cnt);

    if(s.use_cent)
     [cy cx] = ind2sub(size(phi),single(C));
   
     f_cent = ((s.centroid(1) - s.predict_cent(1)) * (cy - s.centroid(1))/s.in.cnt) + ...
             ((s.centroid(2) - s.predict_cent(2)) * (cx - s.centroid(2))/s.in.cnt);
     
     f_y = f_y - (max(abs(f_y))/4) * f_cent/max(abs(f_cent));
     
    end
    
    %%-- (2) texture
%     normalize = @(img) img / (sqrt(dot(img,img))+eps);
%     I = normalize(img(C));
%     I_hat = normalize(s.img_hat(C));
%     inside = find(phi <= 0);
%     in_prod = dot(normalize(img(inside)), normalize(s.img_hat(inside)));
%     f_tx = in_prod*(I.^2 + I_hat.^2)/2 - (I .* I_hat);

    %%-- (3) combined energy
%     f_tx_max = max(abs(f_tx));
%     f_y_max = max(abs(f_y));
%     f = f_y + 0*(f_y_max)*(f_tx/(f_tx_max+eps));
    f = f_y;

    [dx dy] = gradient(single(phi));
    dx = dx(C); dy = dy(C);
    normalize = sqrt(dx.^2 + dy.^2);
    N = [ dx; dy ] ./ (normalize([1 1],:) + eps);
    fN = f([1 1],:) .* N; % x,y components of speed along normal

    mean_fN = mean(fN,2);
    gamma = norm(mean_fN)/10;
    gN = gamma*fN + (1-gamma)*mean_fN(:,ones([1 length(fN)]));
    g = sum( gN .* N, 1);
    S = int8( g*0 );
    return;


    %%-- compute and normalize the gradient and normal --%%
    dx = Dx(phi); dy = Dy(phi);
    dx = single(dx(C)); dy = single(dy(C));
    normalize = sqrt(dx.^2 + dy.^2);
    N = [ dx; dy ] ./ (normalize([1 1],:) + eps);
    
    %%-- compute H1 gradient --%%
    S = int8( sobolev(f, N) );

  end % init_iteration



  function [from, to] = update(from, to, p) % vectorized
    f_unpacked = from.mean * from.cnt;
    t_unpacked = to.mean * to.cnt;
    % add
    to.cnt = to.cnt + numel(p);
    to.mean = (t_unpacked + sum(s.img_P(p)))/to.cnt;
    % remove
    from.cnt = from.cnt - numel(p);
    if from.cnt == 0,  from.mean = 0;
    else               from.mean = (f_unpacked - sum(s.img_P(p)))/from.cnt; end
  end
  function move_in(p) % vectorized
    [s.out s.in] = update(s.out, s.in, p); % update means, counts
    % update centroid
    [x y] = unpack_point(single(p)');
    s.centroid = (s.centroid*(s.in.cnt-numel(p)) + sum([y x]))/s.in.cnt;
  end
  function move_out(p) % vectorized
    %if(s.in.cnt <= 3), return; end % keep from disappearing
    [s.in s.out] = update(s.in, s.out, p); % update means, counts
    % update centroid
    [x y] = unpack_point(single(p'));
    s.centroid = (s.centroid*(s.in.cnt+numel(p)) - sum([y x]))/s.in.cnt;
  end

  function set_img(new_img)
    %s.img_hat = img;
    img = new_img;
  end

  function [x, y] = unpack_point(p) % vectorized
    x = ceil(p / size(img,1));
    y = mod(p, size(img,1));
    y(y == 0) = size(img,1);
  end
  function ind = pack_point(x, y) % vectorized
    ind = (x-1) * size(img,1) + y;
  end
  
  
  function P = translate(M, g, pad)
    dy = g(1); dx = g(2);
    P = pad * ones(size(M), class(pad));
    if dy <= 0
      if dx <= 0
        P(1:end+dy,1:end+dx) = M(1-dy:end,1-dx:end);
      else % 0 < dx
        P(1:end+dy,1+dx:end) = M(1-dy:end,1:end-dx);
      end
    else % 0 < dy
      if dx <= 0
        P(1+dy:end,1:end+dx) = M(1:end-dy,1-dx:end);
      else % 0 < dx
        P(1+dy:end,1+dx:end) = M(1:end-dy,1:end-dx);
      end
    end
  end % translate

end
