function fun(ResultTable, Data, resolution)
  channels = Data.O.Channels;
  cell_size_channel = find(ismember(channels, 'SE'));
  STRADa_channel = find(ismember(channels, 'STRADa'));
  
  cell_size = ResultTable.CInt(:,cell_size_channel);
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);

  STRADa_cyto = STRADa_total - STRADa_nuc;
  STRADa_nuc_scaled_for_size = STRADa_nuc./cell_size;
  STRADa_cyto_scaled_for_size = STRADa_cyto./cell_size;

  %% Heat Map
  figure('Position',[251,    91,   1483,   708]);
  hold on
  
  %% CYTO
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size,1.0),prctile(cell_size,99.0),resolution);
  Y = linspace(prctile(STRADa_cyto_scaled_for_size,1.0),prctile(STRADa_cyto_scaled_for_size,99.0),resolution);
  [x, y] = meshgrid(X,Y);

  % Plot Heat Map
  two_dimensions = [ cell_size  STRADa_cyto_scaled_for_size ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  contour(X,Y,JPD,20,'LineColor','b')
  
  %% NUC
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size,1.0),prctile(cell_size,99.0),resolution);
  Y = linspace(prctile(STRADa_nuc_scaled_for_size,1.0),prctile(STRADa_nuc_scaled_for_size,99.0),resolution);
  [x, y] = meshgrid(X,Y);

  % Plot Heat Map
  two_dimensions = [ cell_size  STRADa_nuc_scaled_for_size ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  contour(X,Y,JPD,10,'LineColor','r')

  
  % Make the plot pretty
  xlabel('Cell Size (SE)');ylabel('STRADa / Cell Size (SE)');
  legend({'Nuclear STRADa','Cytoplasmic STRADa'})
  xlim([prctile(cell_size, 1.0) prctile(cell_size, 99.0)]);
  both_STRA = [STRADa_cyto_scaled_for_size STRADa_nuc_scaled_for_size];
  ylim([prctile(both_STRA(:), 1.0) prctile(both_STRA(:), 99.0)]);

  % Save the plot to disk
  filename = sprintf('contour_STRADa_scaled_for_size.png');
  saveas(gcf,['plots\' filename]);
end

