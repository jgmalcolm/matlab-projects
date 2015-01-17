function fn = weight_tensor(q)
  fn = @weights;

  function w = weights(p_img, p)
    fprintf('E=%f\n', dist_tensor(p,q));

    Z = feature_cov(p_img);
    n = size(Z,2);
    
    w = zeros(size(p_img));
    for x = 1:numel(p_img)
      for i = 1:n
        for j = 1:n
          w(x) = w(x) + Z(x,i)*Z(x,j)*q(i,j);
        end
      end
    end
  end

end
