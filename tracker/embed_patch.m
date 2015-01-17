function P = embed_patch(p, c, sz, pad)
  if ~exist('pad', 'var'), pad = p(1); end

  P = pad * ones(sz, class(pad));
  w = (size(p)-1)/2;
  % defaults for big and small images
  yy_start = 1; YY_start = c(1) - w(1);
  yy_end = size(p,1); YY_end = c(1) + w(1);
  xx_start = 1; XX_start = c(2) - w(2);
  xx_end = size(p,2); XX_end = c(2) + w(2);
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
  % embed
  P(YY_start:YY_end, XX_start:XX_end) = p(yy_start:yy_end, xx_start:xx_end);
end
