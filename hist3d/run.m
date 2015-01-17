function run(img, mask, quanta)
  if ~exist('quanta') quanta = 16; end
  sz = size(img);
  img_ = reshape(img, [], 3);
  [fn h] = hist3d(img_(mask,:), quanta);
  
%   obj{1} = hist3d(s.imgs(:,:
%   for i = 1:length(obj)
%   end

  p = fn(img);
  
  p = reshape(p, sz(1:2));
  clf;colormap gray
  subplot(2,3,1); imagesc(img); axis image off;
  subplot(2,3,2); imagesc(mask); axis image off;
  subplot(2,3,3); imagesc(p); axis image off;

  for i = 1:(256/quanta)
    subplot(8,4,16+i); imagesc(h(:,:,i)); axis image off;
  end
