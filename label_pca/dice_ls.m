function d = dice_ls(lm1, lm2, label_cnt)
% DICE Compute Dice coefficient ([0,1] measure of overlap; higher is better)
  
  p1 = label_probability(lm1, label_cnt);
  p2 = label_probability(lm2, label_cnt);
  
  p1 = reshape(p1, [], label_cnt);
  p2 = reshape(p2, [], label_cnt);

  for i = 1:label_cnt
    A = p1(:,i);
    B = p2(:,i);
    num = sum(min(A, B));
    den = mean([sum(A) sum(B)]);
    d(i) = num / den;
  end
end
