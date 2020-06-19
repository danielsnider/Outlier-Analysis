function fun(ResultTable, Data, resolution)
  channels = Data.O.Channels;
  cell_size_channel = find(ismember(channels, 'SE'));
  STRADa_channel = find(ismember(channels, 'STRADa'));
  cell_size = ResultTable.CInt(:,cell_size_channel);
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;
  STRADa_localization_ratio_nuc_cyto = STRADa_nuc./STRADa_cyto;
  STRADa_localization_ratio_nuc_total = STRADa_nuc./STRADa_total;
  STRADa_localization_ratio_cyto_total = STRADa_cyto./STRADa_total;

  all_STRA = [STRADa_localization_ratio_nuc_cyto ...
              STRADa_localization_ratio_nuc_total ...
              STRADa_localization_ratio_cyto_total]; % used for consistent plot axis limits

  figure('Position',[251,    91,   1583,   708]);
  
  %% STRADa_nuc./STRADa_cyto
  subplot(1,3,1)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size,1),prctile(cell_size,99),resolution);
  Y = linspace(prctile(STRADa_localization_ratio_nuc_cyto,1),prctile(STRADa_localization_ratio_nuc_cyto,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size  STRADa_localization_ratio_nuc_cyto ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size) & isfinite(STRADa_localization_ratio_nuc_cyto);
  fit_line = fit(cell_size(idx),STRADa_localization_ratio_nuc_cyto(idx),'poly1');
  plot(fit_line)
  % Make plot pretty
  set(gca,'ydir','normal')
  title('Nuc / Cyto')
  legend(gca,'off');
  set(h,'markersize',.01);
  xlabel('Cell Size (SE)');ylabel('STRADa_n_u_c/STRADa_c_y_t_o');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(all_STRA(:), 1) prctile(all_STRA(:), 99)]);

  % Fix colors in other subplots
  subplot(1,3,2)
  imagesc(X,Y,zeros(size(JPD)));colormap(gca,jet)
  hold on
  subplot(1,3,3)
  imagesc(X,Y,zeros(size(JPD)));colormap(gca,jet)
  hold on
  
  %% STRADa_nuc./STRADa_total
  subplot(1,3,2)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size,1),prctile(cell_size,99),resolution);
  Y = linspace(prctile(STRADa_localization_ratio_nuc_total,1),prctile(STRADa_localization_ratio_nuc_total,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size  STRADa_localization_ratio_nuc_total ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size) & isfinite(STRADa_localization_ratio_nuc_total);
  fit_line = fit(cell_size(idx),STRADa_localization_ratio_nuc_total(idx),'poly1');
  plot(fit_line)
  % Make plot pretty
  set(gca,'ydir','normal')
  title('Nuc / Total')
  legend(gca,'off');
  set(h,'markersize',.01);
  xlabel('Cell Size (SE)');ylabel('STRADa_n_u_c/STRADa_c_y_t_o');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(all_STRA(:), 1) prctile(all_STRA(:), 99)]);
  
  %% STRADa_cyto./STRADa_total
  subplot(1,3,3)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size,1),prctile(cell_size,99),resolution);
  Y = linspace(prctile(STRADa_localization_ratio_cyto_total,1),prctile(STRADa_localization_ratio_cyto_total,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size  STRADa_localization_ratio_cyto_total ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size) & isfinite(STRADa_localization_ratio_cyto_total);
  fit_line = fit(cell_size(idx),STRADa_localization_ratio_cyto_total(idx),'poly1');
  plot(fit_line)
  % Make plot pretty
  set(gca,'ydir','normal')
  title('Cyto / Total')
  legend(gca,'off');
  set(h,'markersize',.01);
  xlabel('Cell Size (SE)');ylabel('STRADa_n_u_c/STRADa_c_y_t_o');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(all_STRA(:), 1) prctile(all_STRA(:), 99)]);

  % Save the plot to disk
  filename = sprintf('3_JPD.png');
  saveas(gcf,['plots\' filename]);
end

