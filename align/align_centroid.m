function imgs_a = align_centroid(imgs)
% ALIGN_CENTROID Align by simple centroid centering
%
% Example:
%  >> imgs_a = align_centroid(imgs);
  
  % --== Construct Set ==--
  n = size(imgs,3);
  for i = 1:n
    img = imgs(:,:,i);
    % centroid
    [rr cc] = find(img);
    r_cent = mean(rr);
    c_cent = mean(cc);
    % adjustment
    a = round(size(img,2)/2) - c_cent; % x-translation
    b = round(size(img,1)/2) - r_cent; % y-translation
    % transform
    [M H R] = transform(a, b, 1, 0);
    imgs_a(:,:,i) = transform_img(img, M, H, R);
  end


  function [M, H, R] = transform(a, b, h, t)
  % TRANSFORM Standard transform matrices
    M = [ 1 0 a; 0 1 b; 0 0 1 ];
    H = [ h 0 0; 0 h 0; 0 0 1 ];
    R = [ cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1 ];
  end % transform

  function trans_img = transform_img(img, M, H, R)
  % TRANSFORM_IMG Affine transformation of an image.  Fill with black
  % intensity pixels (0x00).  Use bilinear interpolation.
  %
  % TRANS_IMG = TRANSFORM_IMG(IMG, M, H, R)
    
    T = maketform('affine', M'*R*H); % note, M transposed, order rearranged
    u = [ 1 size(img,2) ];
    v = [ 1 size(img,1) ];
    trans_img = imtransform(img, T, 'bilinear', ...
                            'UData', u, ...
                            'VData', v, ...
                            'XData', u, ...
                            'YData', v, ...
                            'FillValues', 0); % black
  end % transform_img
    
end % align
