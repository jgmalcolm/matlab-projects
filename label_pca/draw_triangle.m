function draw_triangle(x_, label_cnt)
  switch label_cnt
   case 2
    plot([0 1], [0 0], 'k+-');
    hold on; p = plot(x_,zeros(size(x_)), 'b'); hold off;
    set(p, 'LineWidth', 2);
    if 3 <= size(x_,2) && mod(size(x_,2), 2) == 1
      hold on; plot(x_((end-1)/2), 0, 'r.'); hold off;
    end
   case 3
    C = label_corners(label_cnt);
    plot(C(1,[1:end 1]), C(2,[1:end 1]), 'k');
    hold on; p = plot(x_(1,:), x_(2,:), 'b'); hold off;
    set(p, 'LineWidth', 2);
    if 3 <= size(x_,2) && mod(size(x_,2), 2) == 1
      hold on; plot(x_(1,(end+1)/2), x_(2,(end+1)/2), 'r.'); hold off;
    end
   otherwise
    error(sprintf('unsupported label count (%d)', label_cnt));
  end
end
