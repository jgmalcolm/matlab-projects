function imagesc3(m)
  n = ceil(sqrt(size(m,3)));
  for i = 1:size(m,3)
    subplot(n,n,i); imagesc(m(:,:,i)); axis image off;
  end
end
