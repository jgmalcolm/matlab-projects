function p = grab_patch(P, c, w, pad)
  if ~exist('pad', 'var'), pad = P(1); end

  p = pad * ones(2*w+1, class(pad));
  % defaults for big and small images
  YY_start = c(1) - w(1); yy_start = 1;
  YY_end = c(1) + w(1); yy_end = size(p,1);
  XX_start = c(2) - w(2); xx_start = 1;
  XX_end = c(2) + w(2); xx_end = size(p,2);
  % adjust
  if c(1) - w(1) < 1
    YY_start = 1;
    yy_start = 2 - (c(1) - w(1));
  end
  if c(2) - w(2) < 1
    XX_start = 1;
    xx_start = 2 - (c(2) - w(2));
  end
  if size(P,1) < c(1) + w(1)
    YY_end = size(P,1);
    yy_end = size(p,1) - (c(1) + w(1) - size(P,1));
  end
  if size(P,2) < c(2) + w(2)
    XX_end = size(P,2);
    xx_end = size(p,2) - (c(2) + w(2) - size(P,2));
  end
  % grab
  p(yy_start:yy_end, xx_start:xx_end) = P(YY_start:YY_end, XX_start:XX_end);
end
