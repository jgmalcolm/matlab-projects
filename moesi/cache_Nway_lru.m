function fn = cache_Nway_lru(ways, total_size, line_size)
  fn = cache(total_size, line_size, ways, @init, @access, @find_victim);

  function s = init(line_cnt)
    s = cell([line_cnt 1]);
  end

  function s = access(s, index, way)
    ind = find(s{index+1} ~= way); % avoid duplicates
    s{index+1} = [s{index+1}(ind) way]; % end is MRU
  end
  
  function way = find_victim(s, index)
    if isempty(s{index+1})
      way = 1; % start
    elseif length(s{index+1}) == ways % capacity?
      way = s{index+1}(1); % LRU
    else
      way = length(s{index+1}) + 1; % add
    end
  end

end
