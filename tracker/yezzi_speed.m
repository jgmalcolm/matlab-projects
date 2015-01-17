function h = yezzi_speed
  %%-- parameters --%%
  nu = 12; % histogram smoothing kernel size
  hist_rate = 0; % rate of refreshing stored histograms [1=reset, 0=keep]
  prior_damp = 4; % edge dampening weight on exponential prior
  
  % required
  h.init_iteration = @init_iteration;
  h.move_in = @move_in;
  h.move_out = @move_out;
  % extensions
  h.init = @init;
  h.postprocess = @postprocess;
  h.set_img = @set_img;
  h.get_state = @get_state;
  
  s = [];  % dummy for lexical scoping
  kern = binomfilt(nu); % for smooth histogram sampling


  function create_prior(phi)
    mask = zeros(size(phi));
    mask(phi<=0) = 0.01;
    mask = bwdist(mask) - 2;
    pindx = find(mask > 0);
    nindx = find(mask <=0 );
    mask(pindx) = exp(-mask(pindx) / prior_damp);
    mask(nindx) = 1;
    s.prior = mask;
  end

  function create_P(phi)
    img_int = uint8(s.img(:) + 1); % quantize
    hist  = wkeep(conv(imhist(img_int(phi <= 0)), kern), 256);
    s.in.hist = hist_rate*hist + (1-hist_rate)*s.in.hist;
    s.in.hist = s.in.hist / sum(s.in.hist);
    hist = wkeep(conv(imhist(img_int(phi >  0)), kern), 256);
    s.out.hist = hist_rate*hist + (1-hist_rate)*s.out.hist;
    s.out.hist = s.out.hist / sum(s.out.hist);
  end

  function postprocess(phi)
    %create_prior(phi);
    s.prior = ones(size(phi));
    create_P(phi);
  end

  function s_ = get_state
    s_ = s;
  end
  
  % called at the beginning of each new image
  function init(phi, C)
    img_int = uint8(s.img); % quantize
    if ~isfield(s, 'in') || ~isfield(s, 'out') || ~isfield(s, 'prior')
      %-- first time initialization
      s.in.hist  = wkeep(conv(imhist(img_int(phi <= 0)), kern), 256);
      s.in.hist = s.in.hist / sum(s.in.hist);
      s.out.hist = wkeep(conv(imhist(img_int(phi >  0)), kern), 256);
      s.out.hist = s.out.hist / sum(s.out.hist);
      keyboard
      postprocess(phi);
    end
    
    %-- probabilistic image
    pin = s.in.hist(img_int+1);
    pout = s.out.hist(img_int+1);
    s.img = 256 * pin ./ (pin + pout + eps) .* s.prior;
    
    %-- inside
    idx = find(phi <= 0);
    s.in.area = numel(idx);
    s.in.mean = mean(s.img(idx));
    % centroid of inside
    [yy xx] = ind2sub(size(phi), idx);
    s.centroid = [mean(yy) mean(xx)];
    
    %-- outside
    idx = find(phi > 0);
    s.out.area = numel(idx);
    s.out.mean = mean(s.img(idx));
  end


  

  % called at the beginning of each new iteration
  function S = init_iteration(phi, C)
    subplot(2,1,2);
    imagesc(s.img); axis image off
    hold on; contour(phi - 0.5, [0 0], 'r'); hold off; drawnow;

    u = s.in.mean;
    v = s.out.mean;
    I = s.img(C);

%     f = (u-v) * ((img(C)-v)/s.out.area + (img(C)-u)/s.in.area); % Yezzi
%     f = f * numel(phi);

    f = (u-v) * ((I-v) + (I-u)); % Chan-Vese

    S = int8( f );
  end



  % called after each iteration to update statistics
  function move_in(p) % vectorized
    [s.out s.in] = update(s.out, s.in, s.img, p); % update means, counts
    % update centroid
    [xx yy] = ind2sub(size(s.img), single(p));
    s.centroid = (s.centroid*(s.in.area-numel(p)) + sum([yy' xx']))/s.in.area;
  end
  function move_out(p) % vectorized
    [s.in s.out] = update(s.in, s.out, s.img, p); % update means, counts
    [xx yy] = ind2sub(size(s.img), single(p));
    s.centroid = (s.centroid*(s.in.area+numel(p)) - sum([yy' xx']))/s.in.area;
  end

  function set_img(img)
    s.img = img;
  end
end


function [from to] = update(from, to, img, p) % vectorized
  f_unpacked = from.mean * from.area;
  t_unpacked = to.mean * to.area;
  % add
  to.area = to.area + numel(p);
  to.mean = (t_unpacked + sum(img(p)))/to.area;
  % remove
  from.area = from.area - numel(p);
  if from.area == 0,  from.mean = 0;   % TODO: try dividing by (area+epsilon)
  else                from.mean = (f_unpacked - sum(img(p)))/from.area; end
end
