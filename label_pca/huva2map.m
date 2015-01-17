function huva2map(huva_fn)
  huva = read_HUVA(huva_fn);
  if any(huva(:) == 0), huva = huva + 1; end

  %% form lookup table
  huva_labels = unique(huva);
  for i = 1:length(huva_labels)
    lookup_table(huva_labels(i)) = i;
  end
  lookup_table = uint8(lookup_table);
  
  %% transform huva to map
  map = lookup_table(huva);
  
  map_fn = strrep(huva_fn, '.img', '.map');
  fid = fopen(map_fn, 'wb');
  fwrite(fid, map, 'uint8');
  fclose(fid);
end
