function test
  n = 10;
  pts = round(n*rand(2,n)); pts(pts==0) = 1;
  C = sub2ind([n n], pts(1,:), pts(2,:));
  C_ = order_curve([n n], C);
  [pts_(1,:) pts_(2,:)] = ind2sub([n n], C_);
  
  plot(pts(1,:), pts(2,:), 'r*', ...
       pts_(1,:), pts_(2,:), 'b');
end
