function fn = weight_histogram(q, param)
  bins = length(q);

  fn = @weights;

  bin = @(x) ceil(x/(256/bins)); % which bin in feature space?
  dirac_lookup = dirac(0:bins+1, param.dirac_eps); % pre-compute

  function w = weights(p_img, p)
    %- compute d/dp of Bhattacharyya similarity metric
    d_dp = sqrt(q./(p+eps)); % [bins 1]
    
    %- compute Dirac [numel(img) bins]
    u = repmat(1:bins, [numel(p_img) 1]);
    B = repmat(bin(p_img(:)), [1 bins]);
    U = dirac_lookup(abs(B - u) + 1);
    
    %- form weights
    w = U * d_dp;
  end
end

function d2 = d2_georgiou(p, q)
  g = sum(log(p ./ q));
  d2 = g^2 - sum(log(p ./ q).^2);
end

function rho = rho_georgiou(p, q)
  g = sum(log(p ./ q));
  rho = sum(p ./ q)/exp(g);
end
