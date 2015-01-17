function lmap = smooth(lmap, window, stddev)
% SMOOTH Smooth label map with Gaussian kernel of size nu
% assumes 2D label map (so lmap is 3D)
  

  k = fspecial('gaussian', window, stddev);
  for i = 1:size(lmap, 3)
    lmap(:,:,i) = imfilter(lmap(:,:,i), k, 'replicate');
  end
end

%   k = fspecial('gaussian', window, stddev);
%   sz = size(lmap);
%   for i = 1:size(lmap, 1)
%     m = reshape(lmap(i,:,:), sz(2:3));
%     lmap(i,:,:) = imfilter(m, k, 'replicate');
%   end
