function draw_colorbar(min,max)
  delta = (max - min)/1000;
  xx = min:delta:max;
  
  w = round(.10 * numel(xx));
  xx = xx(ones(1,w),:);
  
  % border
%   bw = 4;
%   [xx(1:bw,:) xx(end-bw:end,:) xx(:,1:bw) xx(:,end-bw:end)] = deal(0);
  

  imagesc(xx, [min max]); axis off;
%   set(gca, 'XTick', [], ...
%            'YTick', [1 size(xx,1)], ...
%            'YTickLabel', ['max';'min'], ...
%            'FontSize', 20);
%   whitebg;
end
