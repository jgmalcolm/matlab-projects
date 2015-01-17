function box_sz = crop_HUVA
  fs = dir('*.img');
  n = length(fs);
  
  %% determine box containing all nonzero elements
  for i = 1:n
    fprintf('processing %s...\n', fs(i).name);
    vol = read_HUVA(fs(i).name);

    % grab size
    sz(:,i) = size(vol);

    % determine min/max index for non-zero elements along each dimension
    for j = 1:3
      [box_min(j,i) box_max(j,i)] = span(vol, j);
    end
  end
  
  %% determine maximum non-zero span along each dimension
  for i = 1:3
    box_sz(i) = max(box_max(i,:) - box_min(i,:));
  end
  % enlarge by a few pixels
  box_sz = box_sz + 5;
  % ensure even
  ind_odd = find(mod(box_sz, 2));
  box_sz(ind_odd) = box_sz(ind_odd) + 1;
  
  %% extract that size sub-volume and write to disk
  center = round((box_max + box_min)/2);

  %% HACK: retain first two dimensions to keep things simple
  box_sz(1:2) = 256;
  center(1:2,:) = 128;

  for i = 1:n
    fprintf('writing out %s...\n', fs(i).name);
    vol = read_HUVA(fs(i).name);
    sub_vol = extract(vol, center(:,i), box_sz);
    write_HUVA(['cropped/' fs(i).name], sub_vol);
  end

end


function [box_min box_max] = span(v, dim)
  v = permute(v, [dim setxor(dim, 1:3)]);
  %% start
  for i = 1:size(v,1)
    slice = v(i,:,:);
    if any(slice(:) ~= 0)
      box_min = i;
      break;
    end
  end
  %% end
  for i = size(v,1):-1:1
    slice = v(i,:,:);
    if any(slice(:) ~= 0)
      box_max = i;
      break;
    end
  end
end

function v_ = extract(v, c, sz)
  % compute extent (top, bottom, left, right)
  rr = (c(1) - sz(1)/2 + 1):(c(1) + sz(1)/2);
  rr(rr < 1) = 1; rr(size(v,1) < rr) = size(v,1);
  cc = (c(2) - sz(2)/2 + 1):(c(2) + sz(2)/2);
  cc(cc < 1) = 1; cc(size(v,2) < cc) = size(v,2);
  dd = (c(3) - sz(3)/2 + 1):(c(3) + sz(3)/2);
  dd(dd < 1) = 1; dd(size(v,3) < dd) = size(v,3);
  v_ = v(rr, cc, dd);
end

function write_HUVA(fn, vol)
  fid = fopen(fn, 'wb');
  fwrite(fid, vol, 'short');
  fclose(fid);
end
