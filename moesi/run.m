function run(tr)
  if ~exist('tr'), tr = read_traces; end
  for i = 1:size(tr,2)
    cpu{i}.stats.transfer = zeros([1 size(tr,2)]);
    cpu{i}.stats.invalidations = 0;
    %cpu{i}.cache = cache_Nway_lru(1, 16*1024, 16); % 16kB, 16 byte, direct
    cpu{i}.cache = cache_Nway_lru(2, 64*1024, 128); % 64kB, 128 byte, 2-way LRU
  end

  for t = 1:size(tr{1},2)
    bus = {};
    %%-- processor commands
    for i = 1:size(tr,2)
      if ~isempty(tr{i}{t})
        op = tr{i}{t};
        is_ew = @() is_elsewhere(op.addr, i, cpu); % delay execution
        [bus cpu{i}] = cpu_issue(bus, cpu{i}, i, op, is_ew);
      end
    end

    %%-- bus commands
    for b = bus
      for i = setxor(1:size(tr,2), b{1}.id)  % all but my own bus commands
        cpu{i} = bus_issue(cpu{i}, i, b{1});
      end
    end
  end
  
  %%-- statistics
  for i = 1:size(tr,2)
    print_stats(cpu{i}.stats);
    print_cache(cpu{i}.cache);
  end

end


function print_stats(s)
  s
end
function print_cache(c)
  counts = zeros([1 5]);
  function update(l)
    if isempty(l)
      counts(5) = counts(5) + 1;
    else
      switch l.state
       case 'M'; counts(1) = counts(1) + 1;
       case 'O'; counts(2) = counts(2) + 1;
       case 'E'; counts(3) = counts(3) + 1;
       case 'S'; counts(4) = counts(4) + 1;
       case 'I'; counts(5) = counts(5) + 1;
      end
    end
  end
  c.iterate(@update);
  counts
end


function cpu = bus_issue(cpu, id, bus)
  state = cpu.cache.get_state(bus.addr);
  update = @(state) cpu.cache.update(bus.addr, state);
  switch state
   case 'M' % modified
    if strcmp(bus.cmd, 'BusRd')
      update('O');
      cpu.stats.transfer(bus.id) = cpu.stats.transfer(bus.id) + 1;
    elseif strcmp(bus.cmd, 'BusRdX')
      update('I');
      cpu.stats.invalidations = cpu.stats.invalidations + 1;
    end
   case 'O' % owner
    if strcmp(bus.cmd, 'BusRdX')
      update('I');
      cpu.stats.invalidations = cpu.stats.invalidations + 1;
    elseif strcmp(bus.cmd, 'BusUpgr')
      update('I');
      cpu.stats.invalidations = cpu.stats.invalidations + 1;
    elseif strcmp(bus.cmd, 'BusRd')
      update('O');
      cpu.stats.transfer(bus.id) = cpu.stats.transfer(bus.id) + 1;
    end
   case 'E' % exclusive
    if strcmp(bus.cmd, 'BusRd')
      update('S');
    elseif strcmp(bus.cmd, 'BusRdX')
      update('I');
      cpu.stats.invalidations = cpu.stats.invalidations + 1;
    end
   case 'S' % shared
    if strcmp(bus.cmd, 'BusRd')
      update('S');
    elseif strcmp(bus.cmd, 'BusRdX') || strcmp(bus.cmd, 'BusUpgr')
      update('I');
      cpu.stats.invalidations = cpu.stats.invalidations + 1;
    end
   case 'I' % invalid
  end
end

function [bus cpu] = cpu_issue(bus, cpu, id, op, is_elsewhere)
  state = cpu.cache.get_state(op.addr);
  switch state
   case 'M' % modified
   case 'O' % owner
    if op.is_write
      state = 'M';
      bus{end+1} = struct('id', id, 'addr', op.addr, 'cmd', 'BusUpgr');
    end
   case 'E' % exclusive
    if op.is_write
      state = 'M';
    end
   case 'S' % shared
    if op.is_write
      state = 'M';
      bus{end+1} = struct('id', id, 'addr', op.addr, 'cmd', 'BusUpgr');
    end
   case 'I' % invalid
    if op.is_write
      state = 'M';
      bus{end+1} = struct('id', id, 'addr', op.addr, 'cmd', 'BusRdX');
    else
      if is_elsewhere() % delayed execution of check
        state = 'S';
      else
        state = 'E';
      end
      bus{end+1} = struct('id', id, 'addr', op.addr, 'cmd', 'BusRd');
    end
  end
  cpu.cache.update(op.addr, state);
end

function r = is_elsewhere(addr, id, cpu)
  for i = setxor(1:size(cpu,2), id) % skip self
    state = cpu{i}.cache.get_state(addr);
    if state ~= 'I'
      r = 1; % found elsewhere
      return
    end
  end
  r = 0; % default: not found
end
