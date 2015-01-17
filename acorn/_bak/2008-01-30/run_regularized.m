function run_regularized
  %ex = 'SYNTH1'
  %ex = 'ISRAEL_movie14'
  %ex = 'VAN'
  ex = 'SYNTH2'

  paths
  s = load(['images/' ex ], 'imgs', 'win', 'init');
  
  E_unregularized = landscape(s.imgs{1}, s.win, s.init, ...
                              model_tensor(s.win, @kern_none), ...
                              @dist_L2);
  E_regularized   = landscape(s.imgs{1}, s.win, s.init, ...
                              model_tensor(s.win, @kern_Epanechnikov), ...
                              @dist_L2);
  
  clf; colormap gray

  subplot(2,1,1);
  contour(E_unregularized);
  hold on; plot(s.init(2), s.init(1), 'r.'); hold off;
  axis ij square;

  subplot(2,1,2);
  contour(E_regularized);
  hold on; plot(s.init(2), s.init(1), 'r.'); hold off;
  axis ij square;
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
