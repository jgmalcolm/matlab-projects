function run(test_ind)
  if ~exist('test_ind'), test_ind = 1; end
  
  map = @label_map; unmap = @label_unmap;
%   map = @sdf_map; unmap = @sdf_unmap;

  %s = load('images/FOUR');
  s = load('images/COVARY_TOGETHER');
  %s = load('images/COVARY_SEPARATE');

  [r c n] = size(s.imgs);
  label_cnt = numel(unique(s.imgs));
  
  %% training set
  train_set = setxor(test_ind, 1:n);
  %train_set = 1:n;  %% uncomment to include all examples in training set
  

  % map data into appropriate domain
  for i = 1:n
    x_(:,i) = f(s.imgs(:,:,i));
  end

  % test image
  clf; colormap gray; subplot(3,3,1);
  test_img = s.imgs(:,:,test_ind);
  test_img = f_inv(f(test_img));
  imagesc(test_img); title('test image'); axis image off;
  
  % kernel pre-image
  fn = shape_kpca(x_(:,train_set));
  x_test = fn.preimage(x_(:,test_ind));
  img = f_inv(x_test);
  subplot(3,3,3);
  imagesc(img); axis image off;
%   title(sprintf('kernel pre-image (%.1f%% wrong)', ...
%                 misslabeled(test_img, img)*100));
  
  % linear pre-image
  fn = shape_pca(x_(:,train_set));
  x_test = fn.preimage(x_(:,test_ind));
  img = f_inv(x_test);
  subplot(3,3,2);
  imagesc(img); axis image off;
%   title(sprintf('linear pre-image (%.1f%% wrong)', ...
%                 misslabeled(test_img, img)*100));

  %% modes of variation
  scales = -1:.5:1;
%   scales = -2:2;
  i = 1;
  for mode = 1:4
    for scale = scales
      subplot(6,length(scales),2*length(scales)+i);
      if isequal(unmap, @sdf_unmap)
        x_hat = fn.preimage_mode(mode, scale);
        x_hat = reshape(x_hat, [size(img) label_cnt-1]);
        % draw each contour
        imagesc(zeros(size(img)));
        hold on;
        colors = ' rkbmy';
        for lbl = 2:label_cnt
          contour(x_hat(:,:,lbl-1), [0 0], colors(lbl));
        end
        hold off;
      else
        imagesc(img_mode(mode, scale));
      end
      axis image off;
      if scale == 0    title('\mu');
      else             title(sprintf('%.1f\\sigma', scale)); end
      i = i + 1;
    end
  end

  function x_ = f(img)
    x_ = map(img(:), label_cnt);
    x_ = x_(:);
  end
  
  function img = f_inv(x_)
    % which label closest
    img_lbl = unmap(x_, label_cnt);
    img_lbl = reshape(img_lbl, [r c]);
%     img = img_lbl;
%     return
    
    % how confident
    img_p = label_probability(x_, label_cnt);
    img_p = img_p(sub2ind(size(img_p), img_lbl(:), (1:size(img_p,2))'));
    img_p = reshape(img_p, [r c]);
    
    % scale color of label by confidence
    colors = 255*[1 0 0; 0 1 0; 0 0 1];
    img = colors(img_lbl, :);
    img = reshape(img, [r c 3]) .* img_p(:,:,[1 1 1]);
    img = uint8(img);
  end
  
  function p = misslabeled(truth, test)
    p = sum(test(:) ~= truth(:)) / numel(test);
  end
  
  function img = img_mode(mode, scale)
    x_hat = fn.preimage_mode(mode, scale);
    img = f_inv(x_hat);
  end
  
  function d2 = corner_norm(a, b)
    error('unfinished');
  end

  function d2 = norm_L1(a,b)
    d2 = sum(abs(a-b));
  end
  
end
