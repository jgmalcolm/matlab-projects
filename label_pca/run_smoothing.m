function map_ = run_smoothing(map, windows, sigma)
  %% parameters
  p.is_writemode = true;
  p.is_colorized = true;
  %p.map = @scalar_map; p.unmap = @scalar_unmap;
  %p.map = @s_map; p.unmap = @s_unmap;
  %p.map = @label_map; p.unmap = @label_unmap;
  p.map = @tsai_map; p.unmap = @tsai_unmap;

  [r c] = size(map);
  label_cnt = single(max(unique(map)));
  
  for i = 1:length(windows)
    %% map, smooth, and un-map
    map_ = f_inv(f_smooth(f(map)));
    imagesc(map_); axis image off;
    if p.is_writemode
      print('-deps2c', sprintf('figs/s%02d.eps', windows(i)));
    end
  end

  function x_ = f_smooth(x_)
    x_ = reshape(x_, [r c numel(x_)/(r*c)]);
    x_ = smooth(x_, windows(i), sigma);
  end

  function x_ = f(m)
    x_ = p.map(m(:), label_cnt);
    x_ = x_(:);
  end
  
  function m = f_inv(x_)
    % which label closest?
    if p.is_colorized
      m = uint8(m_inv_c(x_, p.unmap, label_cnt));
      m = reshape(m, [r c 3]);
    else
      m = m_inv(x_, p.unmap, label_cnt);
      m = reshape(m, [r c]);
    end
  end
end
