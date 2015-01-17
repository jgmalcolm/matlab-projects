function P = play_point(imgs)

  clf; colormap gray;
  imagesc(imgs{1}); axis image off
  [x(2) x(1)] = ginput(1);
  
  P = zeros(6, numel(imgs));

  last_img = imgs{1};
  for i = 1:numel(imgs)
    img = imgs{i};

    % process
    p = register(last_img, img);
    x = project(x, p);
    
    P(:,i) = p;

    % draw
    imagesc(img); axis image off; title(i);
    hold on; plot(x(2), x(1), 'r.'); hold off;
    drawnow;
    
    % save
    last_img = img;
  end

end

function x = project(x, P)
  a = P(1); b = P(2); c = P(3);
  d = P(4); e = P(5); f = P(6);
  x(2,:) = a*x(2) + b*x(1) + c;
  x(1,:) = d*x(2) + e*x(1) + f;
end
