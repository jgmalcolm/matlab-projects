function [fn h] = hist3d(data, quanta)
  if ~exist('quanta') quanta = 16; end

  fn = @evaluate;

  sz = [1 1 1]*256/quanta;

  ind = data2ind(data);
  h = accumarray(ind, 1, [prod(sz) 1]) / numel(ind);
  h = reshape(h, sz);

  % smooth bins
%   k = [1 3 1]/4;
%   h = convn(h, k, 'same');
%   h = convn(h, k', 'same');
%   h = convn(h, reshape(k, [1 1 length(k)]), 'same');
  
  function p = evaluate(data)
    ind = data2ind(data);
    p = h(ind);
  end

  function ind = data2ind(data)
    if ndims(data) == 3
      data = reshape(data, [size(data,1)*size(data,2) 3]);
    end
    data = double(floor(single(data) / quanta)) + 1;
    ind = sub2ind(sz, data(:,1), data(:,2), data(:,3));
  end
end
