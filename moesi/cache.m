function fn = cache(total_size, line_size, ways, init, access, find_victim)
  fn.get_state = @get_state;
  fn.update = @update;
  fn.iterate = @iterate;

  % initialize
  line_cnt = total_size / (line_size * ways);
  lines = cell([line_cnt ways]);
  s = init(line_cnt);

  get_index = @(addr) mod(floor(addr / line_size), line_cnt);
  get_tag = @(addr) floor(addr / line_size / line_cnt);
  
  function state = get_state(addr)
    state = 'I'; % invalid
    ind = get_index(addr);
    for i = 1:ways
      if ~isempty(lines{ind+1,i}) && lines{ind+1,i}.tag == get_tag(addr)
        state = lines{ind+1,i}.state;
        return;
      end
    end
  end

  function update(addr, state)
    ind = get_index(addr);
    way = NaN;
    for i = 1:ways
      if ~isempty(lines{ind+1,i}) && lines{ind+1,i}.tag == get_tag(addr)
        way = i;
        break;
      end
    end
    if isnan(way)
      way = find_victim(s, ind);
    end
    lines{ind+1,way}.tag = get_tag(addr);
    lines{ind+1,way}.state = state;
    s = access(s, ind, way);
  end
  
  function iterate(fn)
    for i = 1:numel(lines)
      fn(lines{i});
    end
  end
end
