function fun(ResultTable, Data, resolution)
  channels = Data.O.Channels;
  DAPI_channel = find(ismember(channels, 'DAPI'));
  geminin_channel = find(ismember(channels, 'Geminin'));
  DAPI = ResultTable.NInt(:,DAPI_channel);
  geminin = ResultTable.NInt(:,geminin_channel);

  % Find cell stages
  [idxEG1,idxLG1,idxG1S,idxS,idxG2] = FindStages_v2(DAPI,log(geminin),'image');

  % combine LG1 and G1S into one group (LG1) because G1S wasn't detected well enough
  idx_LG1 = idx_LG1 | idx_G1S;
  idxG1S = []; % delete so no mistake is made

end