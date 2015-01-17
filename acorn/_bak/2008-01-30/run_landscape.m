function E = run_landscape
  %ex = 'SYNTH1'
  %ex = 'ISRAEL_movie14'
  %ex = 'VAN'
  ex = 'SYNTH2'

  paths
  s = load(['images/' ex ], 'imgs', 'win', 'init');
  
  span = -15:15;
  img = s.imgs{1};
  s.imgs{1} = img(s.init(1) + span, s.init(2) + span);
  s.init = repmat((length(span)+1)/2, [1 2]);

%   E = landscape(s.imgs{1}, s.win, s.init, ...
%                 model_tensor(s.win, @kern_none), ...
%                 @sim_tensor);
  E   = landscape(s.imgs{1}, s.win, s.init, ...
                  model_histogram(size(s.imgs{1})), ...
                  @sim_histogram_y);
  
  clf; colormap gray
  imagesc(E);
  hold on; plot(s.init(2), s.init(1), 'r.'); hold off;
  axis ij square; title('tensor');
end

function k = kern_none(x)
  k = ones(size(x));
end
function k = kern_Epanechnikov(x) % Epanechnikov (k' => constant)
  d = 2; c_d = pi; % circle in 2D
  k = zeros(size(x));
  ind = find(x < 1);
  k(ind) = (d+2)*(1-x(ind))/(2*c_d);
end
