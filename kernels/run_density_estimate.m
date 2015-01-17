function run_density_estimate
  s = load('images/PLANES', 'masks');
  
  for i = 1:length(s.masks)
    x_train(:,i) = f(s.masks{i});
  end
  
  shape_fn = shape_lpca(x_train);
  n = shape_fn.get_basis_cnt();
  [phi_basis phi_sigma phi_mean] = shape_fn.get_basis();

  alpha = zeros([n 1]);
  alpha(2) = 0.2
  x_ = phi_basis * alpha + phi_mean;
  x_(x_ < 0) = 0;

%   stab = @(x) (x + abs(x))/2;
%   x_ = stab(phi_basis * alpha) + stab(phi_mean);

  contourf(f_(x_)); axis ij square off;

  
  function x = f(m)
    %x = single(m(:)); % mask
    x = reshape(mask2kernel(m), [numel(m) 1]); % kernel
  end
  function m = f_(x)
    m = reshape(x, size(s.masks{1}));
  end
end
