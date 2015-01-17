function m = smooth2(m, window, stddev)
% SMOOTH Smooth label map with Gaussian kernel of size nu assumes 2D spatial
% map (so map is 3D with label dimension)

  k = fspecial('gaussian', window, stddev);
  for i = 1:size(m, 3)
    m(:,:,i) = imfilter(m(:,:,i), k, 'replicate');
  end
end
