function fun(ResultTable, Data, resolution)
  channels = Data.O.Channels;
  cell_size_channel = find(ismember(channels, 'SE'));
  STRADa_channel = find(ismember(channels, 'STRADa'));
  
  cell_size = ResultTable.CArea(:);
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;

  % Eliminate outliers
  lower_prctile = .1;
  higher_prctile = 99.9;
  cell_size = set_outliers_to(cell_size, lower_prctile, higher_prctile, NaN);
  STRADa_nuc = set_outliers_to(STRADa_nuc, lower_prctile, higher_prctile, NaN);
  STRADa_total = set_outliers_to(STRADa_total, lower_prctile, higher_prctile, NaN);
  STRADa_cyto = set_outliers_to(STRADa_cyto, lower_prctile, higher_prctile, NaN);

  STRADa_nuc_scaled_for_size = STRADa_nuc./cell_size;
  STRADa_cyto_scaled_for_size = STRADa_cyto./cell_size;

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
  contour(X,Y,JPD,9,'LineColor','b')
  
  %% NUC
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size,1.0),prctile(cell_size,99.0),resolution);
  Y = linspace(prctile(STRADa_nuc_scaled_for_size,1.0),prctile(STRADa_nuc_scaled_for_size,99.0),resolution);
  [x, y] = meshgrid(X,Y);

  % Plot Heat Map
  two_dimensions = [ cell_size  STRADa_nuc_scaled_for_size ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  contour(X,Y,JPD,8,'LineColor','r')

  
  % Make the plot pretty
  xlabel('Cell Size');ylabel('STRADa / Cell Size');
  legend({'STRADa in Nucleus','STRADa in Cytoplasm'})

  % Save the plot to disk
  filename = sprintf('contour_STRADa_scaled_for_size.png');
  export_fig(['plots\' filename], '-m2')
end

