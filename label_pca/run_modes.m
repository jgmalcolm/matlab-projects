function run_modes(maps)
  % parameters
  %p.map = @label_map; p.unmap = @label_unmap;
  p.map = @sdf_map; p.unmap = @sdf_unmap;

  p.modes = 1;
%   p.scales = -3:3:3;
  p.scales = -2:2:2;
%   p.scales = -1:1;

  sz = size(maps{1});
  label_cnt = numel(unique(maps{1}));
  
  %% map data into appropriate domain
  for i = 1:length(maps), x_(:,i) = f(maps{i}); end

  % PCA to determine modes
  fn = shape_lpca(x_);
  
  % display modes
  sp = 1;
  clf;
  for mode = p.modes
    for scale = p.scales
      subplot(numel(p.scales), numel(p.modes), sp);
      view_map3(map_mode(mode, scale), label_cnt);
      sp = sp + 1;
    end
  end

  function x_ = f_smooth(x_)  % spatial smoothing
    kern = ones(3,3,3)/3^3;
    for j = 1:label_cnt-1
      x_(:,:,:,j) = convn(x_(:,:,:,j), kern, 'same');
    end
  end


  function x_ = f(m)
    x_ = p.map(m, label_cnt);
    %x_ = f_smooth(x_);
    x_ = x_(:);
  end
  
  function m = f_inv(x_)
    x_ = reshape(x_, [sz label_cnt-1]);
    m = p.unmap(x_, label_cnt);
  end
  
  function m = f_inv_probability(x_)
    % how confident
    m_lbl = p.unmap(x_, label_cnt);
    m = label_probability(x_, label_cnt);
    m = m(sub2ind(size(m), m_lbl(:), (1:size(m,2))'));
    m = reshape(m, [r c]);
  end
  
  function m = map_mode(mode, scale)
    x_hat = fn.preimage_mode(mode, scale);
    m = f_inv(x_hat);
    %m = f_inv_probability(x_hat);
  end
end
