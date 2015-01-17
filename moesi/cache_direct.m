function fn = cache_direct(total_size, line_size)
  fn = cache_Nway_lru(1, total_size, line_size);
