function fn = filter_fitts(ref, w)
% FILTER_FITTS Matched normalized cross-correlation filter.
%
% FN = FILTER_FITTS(REF, W) Create filter initializing with reference gate and
% update weight.
%
% Example:
%  >> A = gallery('invhess',20) * gallery('chebvand',20)^2;
%  >> ref = A(10:16, 12:16);
%  >> fn = filter_fitts(ref, 1/8);
%  >> y = fn(A); % [17 10]

  update_ref = filter_ra(w);
  update_ref(ref); % initialize

  function y = filter(u)
    %%-- (1) match filter convolution --%%
    c = conv2(u - mean(u(:)), ref - mean(ref(:)), 'same');
    [val ind] = max(c(:));
    [y(1) y(2)] = ind2sub(size(u), ind);
    %%-- (2) update reference --%%
    p = grab_patch(u, y, (size(ref)-1)/2, mean(ref(:)));
    ref = update_ref(p);
    
    subplot(1,3,2); imagesc(ref); title('Fitts ref');
    subplot(1,3,3); imagesc(u); title('Fitts signal');
    hold on; plot(y(2), y(1), 'r.'); hold off;
  end
  
  fn = @filter; % return closure reference
end
