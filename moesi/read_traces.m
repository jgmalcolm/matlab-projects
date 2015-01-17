function tr = read_traces
  tr{1} = get_trace('traces/p0.tr');
  tr{2} = get_trace('traces/p1.tr');
  tr{3} = get_trace('traces/p2.tr');
  tr{4} = get_trace('traces/p3.tr');

  % pad to end
  t_max = max([length(tr{1}) length(tr{2}) length(tr{3}) length(tr{4})]);
  if length(tr{1}) < t_max, tr{1}{t_max} = []; end
  if length(tr{2}) < t_max, tr{2}{t_max} = []; end
  if length(tr{3}) < t_max, tr{3}{t_max} = []; end
  if length(tr{4}) < t_max, tr{4}{t_max} = []; end
end


function tr = get_trace(filename)
  fid = fopen(filename, 'r');
  data = fscanf(fid, '%d %d 0x%x\n', [3 Inf])';
  fclose(fid);
  
  for i = 1:size(data,1)
    timestamp = data(i,1);
    tr{timestamp}.is_write = data(i,2);
    tr{timestamp}.addr = data(i,3);
  end
end
