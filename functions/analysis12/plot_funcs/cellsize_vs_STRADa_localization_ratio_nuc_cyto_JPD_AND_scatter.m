function fun(ResultTable, Data, resolution, special)
  channels = Data.O.Channels;
  STRADa_channel = find(ismember(channels, 'STRADa'));
  DAPI_channel = find(ismember(channels, 'DAPI'));
  geminin_channel = find(ismember(channels, 'Geminin'));
  DAPI = ResultTable.NInt(:,DAPI_channel);
  geminin = ResultTable.NInt(:,geminin_channel);
  cell_size = ResultTable.CArea(:);
  nuc_size = ResultTable.NArea(:);
  cyto_size = cell_size - nuc_size;
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;
  STRADa_nuc = STRADa_nuc./nuc_size;
  STRADa_total = STRADa_total./cell_size;
  STRADa_cyto = STRADa_cyto./cyto_size;
  STRADa_localization_ratio_nuc_cyto = STRADa_nuc./STRADa_cyto;
  STRADa_localization_ratio_nuc_total = STRADa_nuc./STRADa_total;
  STRADa_localization_ratio_cyto_total = STRADa_cyto./STRADa_total;

  % Eliminate outliers
  lower_prctile = 1;
  higher_prctile = 95;
  cell_size = set_outliers_to(cell_size, lower_prctile, higher_prctile, NaN);
  STRADa_nuc = set_outliers_to(STRADa_nuc, lower_prctile, higher_prctile, NaN);
  STRADa_total = set_outliers_to(STRADa_total, lower_prctile, higher_prctile, NaN);
  STRADa_cyto = set_outliers_to(STRADa_cyto, lower_prctile, higher_prctile, NaN);
  STRADa_localization_ratio_nuc_cyto = set_outliers_to(STRADa_localization_ratio_nuc_cyto, lower_prctile, higher_prctile, NaN);
  STRADa_localization_ratio_nuc_total = set_outliers_to(STRADa_localization_ratio_nuc_total, lower_prctile, higher_prctile, NaN);
  STRADa_localization_ratio_cyto_total = set_outliers_to(STRADa_localization_ratio_cyto_total, lower_prctile, higher_prctile, NaN);

  % Find cell stages
  [idx_EG1,idx_LG1,idx_G1S,idx_S,idx_G2] = FindStages_v2(DAPI,log(geminin));

  % combine LG1 and G1S into one group (LG1) because G1S wasn't detected well enough
  idx_LG1 = idx_LG1 | idx_G1S;
  idx_G1S = []; % delete so no mistake is made


  figure('Position',[251,    91,   883,   708]);
  hold on
  
  %% STRADa_nuc./STRADa_cyto
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size,1),prctile(cell_size,99),resolution);
  Y = linspace(prctile(STRADa_localization_ratio_nuc_cyto,1),prctile(STRADa_localization_ratio_nuc_cyto,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size  STRADa_localization_ratio_nuc_cyto ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  if strcmp(special, 'log')
    JPD = log(JPD);
  end
  imagesc(X,Y,JPD);colormap(gca,jet)
  % Make plot pretty
  set(gca,'ydir','normal')
  title('STRADa Localization Ratio vs Cell Area')
  legend(gca,'off');
  xlabel('Cell Size (area)');ylabel('STRADa_n_u_c/STRADa_c_y_t_o (average intensity)');

  % Choose X and Y data to plot
  x_name = 'cell_size';
  y_name = 'STRADa_localization_ratio_nuc_cyto';

  X=eval(x_name);
  Y=eval(y_name);

  % Subset of data so it is not cluttered
  X=X(1:25:end);
  Y=Y(1:25:end);

  if strcmp(special, 'log')
    Y = log(Y);
  end
  h=scatter(X,Y,5,'filled', 'MarkerFaceAlpha',0.35,'MarkerFaceColor','red')
  % set(h,'markersize',3);

  set(gca,'Color',[0 0 0.5608])

  % Save the plot to disk
  filename = sprintf('plots/cellsize_vs_STRADa_localization_ratio_nuc_cyto_JPD_AND_scatter_%s.png',special);
  export_fig(filename, '-m4')
end

