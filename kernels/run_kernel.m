function run_kernel
  s = load('images/PLANES', 'masks');
  
  for i = 1:length(s.masks)
    x_train(:,i) = f(s.masks{i});
  end
  
  shape_fn = shape_lpca(x_train);
  basis_cnt = shape_fn.get_basis_cnt();
  
  alpha = zeros([basis_cnt 1]);
  x_ = shape_fn.preimage(mean(x_train, 2));
  m = f_(x_);
  m(m < 0) = 0;
  contourf(m); axis ij square off;

  
  function x = f(m)
    %x = double(m(:)); % mask
    x = reshape(mask2kernel(m), [numel(m) 1]); % kernel
  end
  function m = f_(x)
    m = reshape(x, size(s.masks{1}));
  end
end
