function d = dice(m1, m2, label_cnt)
% DICE Compute Dice coefficient ([0,1] measure of overlap; higher is better)

  for i = 1:label_cnt
    A = m1(:) == i;
    B = m2(:) == i;
    num = sum(A & B);
    den = mean([sum(A) sum(B)]);
    d(i) = num / den;
  end
end
