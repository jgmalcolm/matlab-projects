function run_tracker(data)

  clf; paths;
  s = tracker_init(data);

  for i = 1:length(data.imgs);
    s = tracker(s, data, i);
  end
end
