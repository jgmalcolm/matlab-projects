function mask = center_mask(mask)
  
  [r c n] = size(mask);

  if length(which('missile_centroid')) == 0, addpath('~gtg233b/functions'); end
  
  center = [r c]/2;

  for i = 1:n
    m = mask(:,:,i);
    [x y] = missile_centroid(m);
    diff = center - [y x];
    fprintf('before: (%f, %f)   %f\n', y, x, norm(diff));

    [M H R] = transform(diff(1), diff(2), 1, 0);
    m_centered = do_transform(m, M, H, R);

    [x y] = missile_centroid(m_centered);
    diff = center - [y x];
    fprintf('after: (%f, %f)   %f\n\n', y, x, norm(diff));
    mask(:,:,i) = m;
  end

  


  function [M, H, R] = transform(a, b, h, t)
  % TRANSFORM Standard transform matrices
    M = [ 1 0 a; 0 1 b; 0 0 1 ];
    H = [ h 0 0; 0 h 0; 0 0 1 ];
    R = [ cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1 ];
  end % transform

  function trans_img = do_transform(img, M, H, R)
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
  
end
