function psi = p2psi(P)
  label_cnt = size(P, ndims(P));
  C = label_corners(label_cnt);
  
  sz = size(P);
  P = reshape(P, [numel(P)/label_cnt label_cnt]);

  eps = 1e-4;
  if norm(sum(P(1,:)) - 1) > eps, error('P should be normalized'); end

  % weighted combination of vertices
  psi = P(:,1) * C(1,:);
  for j = 2:label_cnt
    psi = psi + P(:,j) * C(j,:);
  end
  
  sz(end) = sz(end) - 1;
  psi = reshape(psi, sz);
end
