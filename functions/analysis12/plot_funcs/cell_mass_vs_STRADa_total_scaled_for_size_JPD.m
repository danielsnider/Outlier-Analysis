function fun(ResultTable, Data)
  channels = Data.O.Channels;
  cell_size_channel = find(ismember(channels, 'SE'));
  STRADa_channel = find(ismember(channels, 'STRADa'));
  cell_size = ResultTable.CInt(:,cell_size_channel);

  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;
  STRADa_scaled_for_size = STRADa_total./cell_size;

  %% Heat Map
  figure('Position',[251,    91,   1483,   708]);
  
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size,1),prctile(cell_size,99),200);
  Y = linspace(prctile(STRADa_scaled_for_size,1),prctile(STRADa_scaled_for_size,99),200);
  [x, y] = meshgrid(X,Y);

  % Plot Heat Map
  two_dimensions = [ cell_size  STRADa_scaled_for_size ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[200 200]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  set(gca,'ydir','normal')
  xlabel('Cell Size (SE)');ylabel('Total STRADa / Cell Size (SE)');

  % Save the plot to disk
  filename = sprintf('JPD_STRADa_scaled_for_size.png');
  saveas(gcf,['plots\' filename]);
end

