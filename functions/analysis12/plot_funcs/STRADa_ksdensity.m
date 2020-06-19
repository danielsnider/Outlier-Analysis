function fun(ResultTable, Data)
  figure('Position',[251,    91,   683,   708]);

  channels = Data.O.Channels;
  cell_size_channel = find(ismember(channels, 'SE'));
  STRADa_channel = find(ismember(channels, 'STRADa'));
  cell_size = ResultTable.CInt(:,cell_size_channel);
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;

  subplot(2,2,1)
  [fx,sx]=ksdensity(cell_size, linspace(prctile(cell_size, .5), prctile(cell_size, 95.5), 1000));
  plot(sx,fx);
  title('Cell Size (SE)')
  subplot(2,2,2)
  [fx,sx]=ksdensity(STRADa_total, linspace(prctile(STRADa_total, .5), prctile(STRADa_total, 95.5), 1000));
  plot(sx,fx);
  title('Total STRADa')
  subplot(2,2,3)
  [fx,sx]=ksdensity(STRADa_nuc, linspace(prctile(STRADa_nuc, .5), prctile(STRADa_nuc, 95.5), 1000));
  plot(sx,fx);
  title('Nuclear STRADa')
  subplot(2,2,4)
  [fx,sx]=ksdensity(STRADa_cyto, linspace(prctile(STRADa_cyto, .5), prctile(STRADa_cyto, 95.5), 1000));
  plot(sx,fx);
  title('Cyto STRADa')

  % Save the plot to disk
  filename = sprintf('STRADa_ksdensity.png');
  saveas(gcf,['plots\' filename]);

end