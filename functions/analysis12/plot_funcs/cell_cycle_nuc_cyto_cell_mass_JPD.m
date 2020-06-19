function fun(ResultTable, Data, resolution)
  channels = Data.O.Channels;
  cell_size_channel = find(ismember(channels, 'SE'));
  STRADa_channel = find(ismember(channels, 'STRADa'));
  DAPI_channel = find(ismember(channels, 'DAPI'));
  geminin_channel = find(ismember(channels, 'Geminin'));
  DAPI = ResultTable.NInt(:,DAPI_channel);
  geminin = ResultTable.NInt(:,geminin_channel);
  cell_size = ResultTable.CInt(:,cell_size_channel);
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;
  STRADa_localization_ratio_nuc_cyto = STRADa_nuc./STRADa_cyto;
  STRADa_localization_ratio_nuc_total = STRADa_nuc./STRADa_total;
  STRADa_localization_ratio_cyto_total = STRADa_cyto./STRADa_total;

  % Find cell stages
  [idx_EG1,idx_LG1,idx_G1S,idx_S,idx_G2] = FindStages_v2(DAPI,log(geminin));

  % combine LG1 and G1S into one group (LG1) because G1S wasn't detected well enough
  idx_LG1 = idx_LG1 | idx_G1S;
  idx_G1S = []; % delete so no mistake is made


  figure('Position',[251,    91,   1483,   408]);


  % Get stats for each cell stage
  cell_size_EG1 = cell_size(idx_EG1);
  cell_size_LG1 = cell_size(idx_LG1);
  cell_size_S = cell_size(idx_S);
  cell_size_G2 = cell_size(idx_G2);
  STRA_ratio_EG1 = STRADa_localization_ratio_nuc_cyto(idx_EG1);
  STRA_ratio_LG1 = STRADa_localization_ratio_nuc_cyto(idx_LG1);
  STRA_ratio_S = STRADa_localization_ratio_nuc_cyto(idx_S);
  STRA_ratio_G2 = STRADa_localization_ratio_nuc_cyto(idx_G2);
  
  %% STRADa_nuc./STRADa_cyto
  subplot(3,4,1)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size_EG1,1),prctile(cell_size_EG1,99),resolution);
  Y = linspace(prctile(STRA_ratio_EG1,1),prctile(STRA_ratio_EG1,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size_EG1  STRA_ratio_EG1 ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size_EG1) & isfinite(STRA_ratio_EG1);
  fit_line = fit(cell_size_EG1(idx),STRA_ratio_EG1(idx),'poly1');
  plot(fit_line)
  [r p] = corr(cell_size_EG1(idx),STRA_ratio_EG1(idx));
  % Make plot pretty
  set(gca,'ydir','normal')
  subplot_title = sprintf('p=%.2f, r=%.2f', p, r);
  title({'Early G1',subplot_title})
  legend(gca,'off');
  xlabel('');ylabel({'STRADa','Nuclear / Cytoplasm'});
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio_nuc_cyto(:), 1) prctile(STRADa_localization_ratio_nuc_cyto(:), 99)]);

  %% STRADa_nuc./STRADa_total
  subplot(3,4,2)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size_LG1,1),prctile(cell_size_LG1,99),resolution);
  Y = linspace(prctile(STRA_ratio_LG1,1),prctile(STRA_ratio_LG1,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size_LG1  STRA_ratio_LG1 ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size_LG1) & isfinite(STRA_ratio_LG1);
  fit_line = fit(cell_size_LG1(idx),STRA_ratio_LG1(idx),'poly1');
  plot(fit_line)
  [r p] = corr(cell_size_LG1(idx),STRA_ratio_LG1(idx));
  % Make plot pretty
  set(gca,'ydir','normal')
  subplot_title = sprintf('p=%.2f, r=%.2f', p, r);
  title({'Late G1',subplot_title})
  legend(gca,'off');
  xlabel('');ylabel('');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio_nuc_cyto(:), 1) prctile(STRADa_localization_ratio_nuc_cyto(:), 99)]);
  
  %% STRADa_cyto./STRADa_total
  subplot(3,4,3)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size_S,1),prctile(cell_size_S,99),resolution);
  Y = linspace(prctile(STRA_ratio_S,1),prctile(STRA_ratio_S,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size_S  STRA_ratio_S ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size_S) & isfinite(STRA_ratio_S);
  fit_line = fit(cell_size_S(idx),STRA_ratio_S(idx),'poly1');
  plot(fit_line)
  [r p] = corr(cell_size_S(idx),STRA_ratio_S(idx));
  % Make plot pretty
  set(gca,'ydir','normal')
  subplot_title = sprintf('p=%.2f, r=%.2f', p, r);
  title({'S',subplot_title})
  legend(gca,'off');
  xlabel('');ylabel('');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio_nuc_cyto(:), 1) prctile(STRADa_localization_ratio_nuc_cyto(:), 99)]);

  %% STRADa_cyto./STRADa_total
  subplot(3,4,4)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size_G2,1),prctile(cell_size_G2,99),resolution);
  Y = linspace(prctile(STRA_ratio_G2,1),prctile(STRA_ratio_G2,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size_G2  STRA_ratio_G2 ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size_G2) & isfinite(STRA_ratio_G2);
  fit_line = fit(cell_size_G2(idx),STRA_ratio_G2(idx),'poly1');
  plot(fit_line)
  [r p] = corr(cell_size_G2(idx),STRA_ratio_G2(idx));
  % Make plot pretty
  set(gca,'ydir','normal')
  subplot_title = sprintf('p=%.2f, r=%.2f', p, r);
  title({'G2',subplot_title})
  legend(gca,'off');
  xlabel('');ylabel('');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio_nuc_cyto(:), 1) prctile(STRADa_localization_ratio_nuc_cyto(:), 99)]);








  % Get stats for each cell stage
  cell_size_EG1 = cell_size(idx_EG1);
  cell_size_LG1 = cell_size(idx_LG1);
  cell_size_S = cell_size(idx_S);
  cell_size_G2 = cell_size(idx_G2);
  STRA_ratio_EG1 = STRADa_localization_ratio_nuc_total(idx_EG1);
  STRA_ratio_LG1 = STRADa_localization_ratio_nuc_total(idx_LG1);
  STRA_ratio_S = STRADa_localization_ratio_nuc_total(idx_S);
  STRA_ratio_G2 = STRADa_localization_ratio_nuc_total(idx_G2);
  
  %% STRADa_nuc./STRADa_cyto
  subplot(3,4,5)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size_EG1,1),prctile(cell_size_EG1,99),resolution);
  Y = linspace(prctile(STRA_ratio_EG1,1),prctile(STRA_ratio_EG1,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size_EG1  STRA_ratio_EG1 ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size_EG1) & isfinite(STRA_ratio_EG1);
  fit_line = fit(cell_size_EG1(idx),STRA_ratio_EG1(idx),'poly1');
  plot(fit_line)
  [r p] = corr(cell_size_EG1(idx),STRA_ratio_EG1(idx));
  % Make plot pretty
  set(gca,'ydir','normal')
  subplot_title = sprintf('p=%.2f, r=%.2f', p, r);
  title(subplot_title)
  legend(gca,'off');
  xlabel('');ylabel({'STRADa','Nuclear / Total'});
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio_nuc_total(:), 1) prctile(STRADa_localization_ratio_nuc_total(:), 99)]);

  %% STRADa_nuc./STRADa_total
  subplot(3,4,6)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size_LG1,1),prctile(cell_size_LG1,99),resolution);
  Y = linspace(prctile(STRA_ratio_LG1,1),prctile(STRA_ratio_LG1,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size_LG1  STRA_ratio_LG1 ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size_LG1) & isfinite(STRA_ratio_LG1);
  fit_line = fit(cell_size_LG1(idx),STRA_ratio_LG1(idx),'poly1');
  plot(fit_line)
  [r p] = corr(cell_size_LG1(idx),STRA_ratio_LG1(idx));
  % Make plot pretty
  set(gca,'ydir','normal')
  subplot_title = sprintf('p=%.2f, r=%.2f', p, r);
  title(subplot_title)
  legend(gca,'off');
  xlabel('');ylabel('');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio_nuc_total(:), 1) prctile(STRADa_localization_ratio_nuc_total(:), 99)]);
  
  %% STRADa_cyto./STRADa_total
  subplot(3,4,7)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size_S,1),prctile(cell_size_S,99),resolution);
  Y = linspace(prctile(STRA_ratio_S,1),prctile(STRA_ratio_S,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size_S  STRA_ratio_S ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size_S) & isfinite(STRA_ratio_S);
  fit_line = fit(cell_size_S(idx),STRA_ratio_S(idx),'poly1');
  plot(fit_line)
  [r p] = corr(cell_size_S(idx),STRA_ratio_S(idx));
  % Make plot pretty
  set(gca,'ydir','normal')
  subplot_title = sprintf('p=%.2f, r=%.2f', p, r);
  title(subplot_title)
  legend(gca,'off');
  xlabel('');ylabel('');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio_nuc_total(:), 1) prctile(STRADa_localization_ratio_nuc_total(:), 99)]);

  %% STRADa_cyto./STRADa_total
  subplot(3,4,8)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size_G2,1),prctile(cell_size_G2,99),resolution);
  Y = linspace(prctile(STRA_ratio_G2,1),prctile(STRA_ratio_G2,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size_G2  STRA_ratio_G2 ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size_G2) & isfinite(STRA_ratio_G2);
  fit_line = fit(cell_size_G2(idx),STRA_ratio_G2(idx),'poly1');
  plot(fit_line)
  [r p] = corr(cell_size_G2(idx),STRA_ratio_G2(idx));
  % Make plot pretty
  set(gca,'ydir','normal')
  subplot_title = sprintf('p=%.2f, r=%.2f', p, r);
  title(subplot_title)
  legend(gca,'off');
  xlabel('');ylabel('');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio_nuc_total(:), 1) prctile(STRADa_localization_ratio_nuc_total(:), 99)]);








  % Get stats for each cell stage
  cell_size_EG1 = cell_size(idx_EG1);
  cell_size_LG1 = cell_size(idx_LG1);
  cell_size_S = cell_size(idx_S);
  cell_size_G2 = cell_size(idx_G2);
  STRA_ratio_EG1 = STRADa_localization_ratio_cyto_total(idx_EG1);
  STRA_ratio_LG1 = STRADa_localization_ratio_cyto_total(idx_LG1);
  STRA_ratio_S = STRADa_localization_ratio_cyto_total(idx_S);
  STRA_ratio_G2 = STRADa_localization_ratio_cyto_total(idx_G2);
  

  %% STRADa_nuc./STRADa_cyto
  subplot(3,4,9)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size_EG1,1),prctile(cell_size_EG1,99),resolution);
  Y = linspace(prctile(STRA_ratio_EG1,1),prctile(STRA_ratio_EG1,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size_EG1  STRA_ratio_EG1 ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size_EG1) & isfinite(STRA_ratio_EG1);
  fit_line = fit(cell_size_EG1(idx),STRA_ratio_EG1(idx),'poly1');
  plot(fit_line)
  [r p] = corr(cell_size_EG1(idx),STRA_ratio_EG1(idx));
  % Make plot pretty
  set(gca,'ydir','normal')
  subplot_title = sprintf('p=%.2f, r=%.2f', p, r);
  title(subplot_title)
  legend(gca,'off');
  xlabel('Cell Size (SE)');ylabel({'STRADa','Cytoplasm / Total'});
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio_cyto_total(:), 1) prctile(STRADa_localization_ratio_cyto_total(:), 99)]);

  %% STRADa_nuc./STRADa_total
  subplot(3,4,10)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size_LG1,1),prctile(cell_size_LG1,99),resolution);
  Y = linspace(prctile(STRA_ratio_LG1,1),prctile(STRA_ratio_LG1,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size_LG1  STRA_ratio_LG1 ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size_LG1) & isfinite(STRA_ratio_LG1);
  fit_line = fit(cell_size_LG1(idx),STRA_ratio_LG1(idx),'poly1');
  plot(fit_line)
  [r p] = corr(cell_size_LG1(idx),STRA_ratio_LG1(idx));
  % Make plot pretty
  set(gca,'ydir','normal')
  subplot_title = sprintf('p=%.2f, r=%.2f', p, r);
  title(subplot_title)
  legend(gca,'off');
  xlabel('Cell Size (SE)');ylabel('');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio_cyto_total(:), 1) prctile(STRADa_localization_ratio_cyto_total(:), 99)]);
  
  %% STRADa_cyto./STRADa_total
  subplot(3,4,11)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size_S,1),prctile(cell_size_S,99),resolution);
  Y = linspace(prctile(STRA_ratio_S,1),prctile(STRA_ratio_S,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size_S  STRA_ratio_S ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size_S) & isfinite(STRA_ratio_S);
  fit_line = fit(cell_size_S(idx),STRA_ratio_S(idx),'poly1');
  plot(fit_line)
  [r p] = corr(cell_size_S(idx),STRA_ratio_S(idx));
  % Make plot pretty
  set(gca,'ydir','normal')
  subplot_title = sprintf('p=%.2f, r=%.2f', p, r);
  title(subplot_title)
  legend(gca,'off');
  xlabel('Cell Size (SE)');ylabel('');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio_cyto_total(:), 1) prctile(STRADa_localization_ratio_cyto_total(:), 99)]);

  %% STRADa_cyto./STRADa_total
  subplot(3,4,12)
  % Make a meshgrid for the heatmap
  X = linspace(prctile(cell_size_G2,1),prctile(cell_size_G2,99),resolution);
  Y = linspace(prctile(STRA_ratio_G2,1),prctile(STRA_ratio_G2,99),resolution);
  [x, y] = meshgrid(X,Y);
  % Plot Heat Map
  two_dimensions = [ cell_size_G2  STRA_ratio_G2 ];
  [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  imagesc(X,Y,JPD);colormap(gca,jet)
  hold on
  % Fit line
  idx = isfinite(cell_size_G2) & isfinite(STRA_ratio_G2);
  fit_line = fit(cell_size_G2(idx),STRA_ratio_G2(idx),'poly1');
  plot(fit_line)
  [r p] = corr(cell_size_G2(idx),STRA_ratio_G2(idx));
  % Make plot pretty
  set(gca,'ydir','normal')
  subplot_title = sprintf('p=%.2f, r=%.2f', p, r);
  title(subplot_title)
  legend(gca,'off');
  xlabel('Cell Size (SE)');ylabel('');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(STRADa_localization_ratio_cyto_total(:), 1) prctile(STRADa_localization_ratio_cyto_total(:), 99)]);





  % Fix background colors
  set(subplot(3,4,1),'Color',[0 0 0.5608])
  set(subplot(3,4,2),'Color',[0 0 0.5608])
  set(subplot(3,4,3),'Color',[0 0 0.5608])
  set(subplot(3,4,4),'Color',[0 0 0.5608])

  set(subplot(3,4,5),'Color',[0 0 0.5608])
  set(subplot(3,4,6),'Color',[0 0 0.5608])
  set(subplot(3,4,7),'Color',[0 0 0.5608])
  set(subplot(3,4,8),'Color',[0 0 0.5608])

  set(subplot(3,4,9),'Color',[0 0 0.5608])
  set(subplot(3,4,10),'Color',[0 0 0.5608])
  set(subplot(3,4,11),'Color',[0 0 0.5608])
  set(subplot(3,4,12),'Color',[0 0 0.5608])

  % Save the plot to disk
  filename = sprintf('cell_cycle_JPD_nuc_cyto.png');
  saveas(gcf,['plots\' filename]);
end

