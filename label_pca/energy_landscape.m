function [E_tsai E_label] = energy_landscape(m_ref, m_free)
  paths;
  
  p.is_batch = true;
  p.range = -3:3;
  p.energy = @energy_tsai;
  
  label_cnt = numel(unique([m_ref m_free]));
  % our representation
  l_ref = reshape(label_map(m_ref, label_cnt)', [size(m_ref) label_cnt-1]);
  l_free = reshape(label_map(m_free, label_cnt)', [size(m_free) label_cnt-1]);
  % Andy's representation
  t_ref = reshape(tsai_map(m_ref, label_cnt)', [size(m_ref) label_cnt-1]);
  t_free = reshape(tsai_map(m_free, label_cnt)', [size(m_free) label_cnt-1]);
  
  [E_tsai E_label] = deal([]);
  for x = p.range
    for y = p.range
      m = move(t_free, x, y);
      E_tsai(end+1) = p.energy(t_ref, m);
      m = move(l_free, x, y);
      E_label(end+1) = p.energy(l_ref, m);
    end
  end
  
  E_tsai = reshape(E_tsai, [numel(p.range) numel(p.range)]);
  E_label = reshape(E_label, [numel(p.range) numel(p.range)]);
  
  clf;
  if p.is_batch
    plot_landscape(E_tsai);  print('-deps2c', '-r90', 'figs/e_tsai.eps');
    plot_landscape(E_label); print('-deps2c', '-r90', 'figs/e_label.eps');
  else
    subplot(2,1,1); plot_landscape(E_tsai);
    subplot(2,1,2); plot_landscape(E_label);
  end
  

  function plot_landscape(E)
    contourf(-E); axis square
    set(gca, 'XTickLabel', p.range, 'YTickLabel', p.range, ...
             'FontSize', 14);
    [yy xx] = find(E == min(E(:)));
    hold on; h = plot(xx, yy, 'r.'); hold off;
    set(h, 'MarkerSize', 50)
  end
end

function E = energy_tsai(a, b)
  % (I_i - I_j)^2 } / { (I_i + I_j)^2
  E = sum( (a(:) - b(:)).^2 ) / sum( (a(:) + b(:)).^2 );
end

function E = energy_mean(a, b)
  % (I_i - mu)^2
  mu = (a(:) + b(:))/2;
  E = sum( (a(:) - mu).^2 + (b(:) - mu).^2 );
end

function m = move(m, dx, dy)
  if dx < 0
    m = move(shiftL(m), dx+1, dy);
  elseif dx > 0
    m = move(shiftR(m), dx-1, dy);
  elseif dy < 0
    m = move(shiftU(m), dx, dy+1);
  elseif dy > 0
    m = move(shiftD(m), dx, dy-1);
  end
end
