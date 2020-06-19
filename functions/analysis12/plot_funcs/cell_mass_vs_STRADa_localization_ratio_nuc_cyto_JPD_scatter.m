function fun(ResultTable, Data, resolution)
  channels = Data.O.Channels;
  cell_size_channel = find(ismember(channels, 'SE'));
  STRADa_channel = find(ismember(channels, 'STRADa'));
  cell_size = ResultTable.CInt(:,cell_size_channel);
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;
  STRADa_localization_ratio = STRADa_nuc./STRADa_cyto;

  % Scatter plot
  figure('Position',[251,    91,   683,   708]);
  idx = isfinite(cell_size) & isfinite(STRADa_localization_ratio);
  fit_line = fit(cell_size(idx),STRADa_localization_ratio(idx),'poly1');
  h=plot(fit_line,cell_size, STRADa_localization_ratio, '.b');

  % Make plot pretty
  legend(gca,'off');
  set(h,'markersize',.01);
  xlabel('Cell Size (SE)');ylabel('STRADa_n_u_c/STRADa_c_y_t_o');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio, 1) prctile(STRADa_localization_ratio, 99)]);

  % Save the plot to disk
  filename = sprintf('scatter.png',);
  saveas(gcf,['plots\' filename]);

  %% Heat Map
  figure('Position',[251,    91,   683,   708]);
  hold on
  
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size,1),prctile(cell_size,99),resolution);
  Y = linspace(prctile(STRADa_localization_ratio,1),prctile(STRADa_localization_ratio,99),resolution);
  [x, y] = meshgrid(X,Y);

  % Plot Heat Map
  two_dimensions = [ cell_size  STRADa_localization_ratio ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)

  % Make plot pretty
  set(gca,'ydir','normal')
  xlabel('Cell Size (SE)');ylabel('STRADa_n_u_c/STRADa_c_y_t_o');
  legend(gca,'off');

  % Save the plot to disk
  filename = sprintf('JPD.png',PlateMapRow.plate, PlateMapRow.column, char(PlateMapRow.stain2_name), char(exp_row_order_array(exp_num))); % example result p1c8_p21_pS6_JPD.png
  saveas(gcf,['plots\' filename]);
end
