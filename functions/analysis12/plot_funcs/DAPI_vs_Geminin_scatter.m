function fun(ResultTable, Data, resolution)
  channels = Data.O.Channels;
  DAPI_channel = find(ismember(channels, 'DAPI'));
  geminin_channel = find(ismember(channels, 'Geminin'));
  DAPI = ResultTable.NInt(:,DAPI_channel);
  geminin = ResultTable.NInt(:,geminin_channel);

  % Scatter plot
  figure('Position',[251,    91,   1583,   708]);
  scatter(DAPI,log(geminin),'.')
  xlabel('DAPI');ylabel('Geminin');
  xlim([prctile(DAPI, 1) prctile(DAPI, 99)]);
  ylim([prctile(log(geminin), 1) prctile(log(geminin), 99)]);

  % Save the plot to disk
  filename = sprintf('scatter_DAPI_vs_Geminin.png');
  saveas(gcf,['plots\' filename]);
end

