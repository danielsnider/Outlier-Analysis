 function fun(ResultTable, Data)
  channels = Data.O.Channels;
  cell_stages = {'EG1', 'LG1', 'S', 'G2'};
  STRADa_channel = find(ismember(channels, 'STRADa'));
  cell_mass_channel = find(ismember(channels, 'SE'));
  DAPI_channel = find(ismember(channels, 'DAPI'));
  geminin_channel = find(ismember(channels, 'Geminin'));
  EG1_color = [0/255 255/255 150/255]*.9;
  LG1_color = [232/255 227/255 12/255]*.9;
  S_color = [255/255 102/255 0/255]*.9;
  G2_color = [207/255 12/255 232/255]*.9;
  DAPI = ResultTable.NInt(:,DAPI_channel);
  geminin = ResultTable.NInt(:,geminin_channel);
  cell_mass = ResultTable.CInt(:,cell_mass_channel);
  nuc_mass = ResultTable.NInt(:,cell_mass_channel);
  cyto_mass = cell_mass - nuc_mass;
  cell_size = ResultTable.CArea(:);
  nuc_size = ResultTable.NArea(:);
  cyto_size = cell_size - nuc_size;
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;
  STRADa_nuc = STRADa_nuc./nuc_mass;
  STRADa_total = STRADa_total./cell_mass;
  STRADa_cyto = STRADa_cyto./cyto_mass;
  STRADa_localization_ratio_nuc_cyto = STRADa_nuc./STRADa_cyto;
  STRADa_localization_ratio_nuc_total = STRADa_nuc./STRADa_total;
  STRADa_localization_ratio_cyto_total = STRADa_cyto./STRADa_total;

  % % Eliminate outliers
  lower_prctile = 0.1;
  higher_prctile = 99.9;
  cell_mass = set_outliers_to(cell_mass, lower_prctile, higher_prctile, NaN);
  nuc_mass = set_outliers_to(nuc_mass, lower_prctile, higher_prctile, NaN);
  cyto_mass = set_outliers_to(cyto_mass, lower_prctile, higher_prctile, NaN);
  cell_size = set_outliers_to(cell_size, lower_prctile, higher_prctile, NaN);
  nuc_size = set_outliers_to(nuc_size, lower_prctile, higher_prctile, NaN);
  cyto_size = set_outliers_to(cyto_size, lower_prctile, higher_prctile, NaN);
  STRADa_nuc = set_outliers_to(STRADa_nuc, lower_prctile, higher_prctile, NaN);
  STRADa_total = set_outliers_to(STRADa_total, lower_prctile, higher_prctile, NaN);
  STRADa_cyto = set_outliers_to(STRADa_cyto, lower_prctile, higher_prctile, NaN);
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


  % Get stats for each cell stage
  cell_size_EG1 = cell_size(idx_EG1);
  cell_size_LG1 = cell_size(idx_LG1);
  cell_size_S = cell_size(idx_S);
  cell_size_G2 = cell_size(idx_G2);
  cell_mass_EG1 = cell_mass(idx_EG1);
  cell_mass_LG1 = cell_mass(idx_LG1);
  cell_mass_S = cell_mass(idx_S);
  cell_mass_G2 = cell_mass(idx_G2);
  STRADa_total_EG1 = STRADa_total(idx_EG1);
  STRADa_total_LG1 = STRADa_total(idx_LG1);
  STRADa_total_S = STRADa_total(idx_S);
  STRADa_total_G2 = STRADa_total(idx_G2);
  STRADa_nuc_EG1 = STRADa_nuc(idx_EG1);
  STRADa_nuc_LG1 = STRADa_nuc(idx_LG1);
  STRADa_nuc_S = STRADa_nuc(idx_S);
  STRADa_nuc_G2 = STRADa_nuc(idx_G2);
  STRADa_cyto_EG1 = STRADa_cyto(idx_EG1);
  STRADa_cyto_LG1 = STRADa_cyto(idx_LG1);
  STRADa_cyto_S = STRADa_cyto(idx_S);
  STRADa_cyto_G2 = STRADa_cyto(idx_G2);
  STRADa_ratio_EG1 = STRADa_localization_ratio_nuc_cyto(idx_EG1);
  STRADa_ratio_LG1 = STRADa_localization_ratio_nuc_cyto(idx_LG1);
  STRADa_ratio_S = STRADa_localization_ratio_nuc_cyto(idx_S);
  STRADa_ratio_G2 = STRADa_localization_ratio_nuc_cyto(idx_G2);


  STRADa_nuc_EG1_mean = nanmean(STRADa_nuc_EG1);
  STRADa_nuc_LG1_mean = nanmean(STRADa_nuc_LG1);
  STRADa_nuc_S_mean = nanmean(STRADa_nuc_S);
  STRADa_nuc_G2_mean = nanmean(STRADa_nuc_G2);
  STRADa_ratio_EG1_mean = nanmean(STRADa_ratio_EG1);
  STRADa_ratio_LG1_mean = nanmean(STRADa_ratio_LG1);
  STRADa_ratio_S_mean = nanmean(STRADa_ratio_S);
  STRADa_ratio_G2_mean = nanmean(STRADa_ratio_G2);
  cell_mass_EG1_mean = nanmean(cell_mass_EG1);
  cell_mass_LG1_mean = nanmean(cell_mass_LG1);
  cell_mass_S_mean = nanmean(cell_mass_S);
  cell_mass_G2_mean = nanmean(cell_mass_G2);



  %% SE vs Ratio
  figure('Position',[680   475   784   623]);
  hold on
  idx = isfinite(cell_mass_EG1) & isfinite(STRADa_ratio_EG1);
  h=plot(cell_mass_EG1(idx),STRADa_ratio_EG1(idx),'.','color',EG1_color);
  idx = isfinite(cell_mass_LG1) & isfinite(STRADa_ratio_LG1);
  h=plot(cell_mass_LG1(idx),STRADa_ratio_LG1(idx),'.','color',LG1_color);
  idx = isfinite(cell_mass_S) & isfinite(STRADa_ratio_S);
  h=plot(cell_mass_S(idx),STRADa_ratio_S(idx),'.','color',S_color);
  idx = isfinite(cell_mass_G2) & isfinite(STRADa_ratio_G2);
  h=plot(cell_mass_G2(idx),STRADa_ratio_G2(idx),'.','color',G2_color);
  idx = isfinite(cell_mass_EG1) & isfinite(STRADa_ratio_EG1);
  fit_line = fit(cell_mass_EG1(idx),STRADa_ratio_EG1(idx),'poly1');
  h=plot(fit_line);
  h.Color = EG1_color .*.7;
  h.LineWidth = 3;
  idx = isfinite(cell_mass_LG1) & isfinite(STRADa_ratio_LG1);
  fit_line = fit(cell_mass_LG1(idx),STRADa_ratio_LG1(idx),'poly1');
  h=plot(fit_line);
  h.Color = LG1_color .*.7;
  h.LineWidth = 3;
  idx = isfinite(cell_mass_S) & isfinite(STRADa_ratio_S);
  fit_line = fit(cell_mass_S(idx),STRADa_ratio_S(idx),'poly1');
  h=plot(fit_line);
  h.Color = S_color .*.7;
  h.LineWidth = 3;
  idx = isfinite(cell_mass_G2) & isfinite(STRADa_ratio_G2);
  fit_line = fit(cell_mass_G2(idx),STRADa_ratio_G2(idx),'poly1');
  h=plot(fit_line);
  h.Color = G2_color .*.7;
  h.LineWidth = 3;
  % Scatter plot by cell cycle MEAN
  h=scatter(cell_mass_EG1_mean,STRADa_ratio_EG1_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',EG1_color);
  h=scatter(cell_mass_LG1_mean,STRADa_ratio_LG1_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',LG1_color);
  h=scatter(cell_mass_S_mean,STRADa_ratio_S_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',S_color);
  h=scatter(cell_mass_G2_mean,STRADa_ratio_G2_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',G2_color);
  h=scatter(cell_mass_G2_mean,STRADa_ratio_G2_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',G2_color);

  legend({'EG1 Cell','LG1 Cell','S Cell','G2 Cell','EG1 Trend','LG1 Trend','S Trend','G2 Trend','EG1 Mean', 'LG1 Mean', 'S Mean', 'G2 Mean'});
  title('Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)')
  xlabel('Cell Mass (SE)');
  ylabel('<- More STRADa in Cytoplasm         Ratio         More STRADa in Nucleus ->');
  axis tight

  filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5 - binned plots/mass-vs-ratio.png');
  export_fig(filename, '-m2 -transparent -nocrop')



  %% SE vs STRADa Nuc
  figure('Position',[680   475   784   623]);
  hold on
  idx = isfinite(cell_mass_EG1) & isfinite(STRADa_nuc_EG1);
  h=plot(cell_mass_EG1(idx),STRADa_nuc_EG1(idx),'.','color',EG1_color);
  idx = isfinite(cell_mass_LG1) & isfinite(STRADa_nuc_LG1);
  h=plot(cell_mass_LG1(idx),STRADa_nuc_LG1(idx),'.','color',LG1_color);
  idx = isfinite(cell_mass_S) & isfinite(STRADa_nuc_S);
  h=plot(cell_mass_S(idx),STRADa_nuc_S(idx),'.','color',S_color);
  idx = isfinite(cell_mass_G2) & isfinite(STRADa_nuc_G2);
  h=plot(cell_mass_G2(idx),STRADa_nuc_G2(idx),'.','color',G2_color);
  idx = isfinite(cell_mass_EG1) & isfinite(STRADa_nuc_EG1);
  fit_line = fit(cell_mass_EG1(idx),STRADa_nuc_EG1(idx),'poly1');
  h=plot(fit_line);
  h.Color = EG1_color .*.7;
  h.LineWidth = 3;
  idx = isfinite(cell_mass_LG1) & isfinite(STRADa_nuc_LG1);
  fit_line = fit(cell_mass_LG1(idx),STRADa_nuc_LG1(idx),'poly1');
  h=plot(fit_line);
  h.Color = LG1_color .*.7;
  h.LineWidth = 3;
  idx = isfinite(cell_mass_S) & isfinite(STRADa_nuc_S);
  fit_line = fit(cell_mass_S(idx),STRADa_nuc_S(idx),'poly1');
  h=plot(fit_line);
  h.Color = S_color .*.7;
  h.LineWidth = 3;
  idx = isfinite(cell_mass_G2) & isfinite(STRADa_nuc_G2);
  fit_line = fit(cell_mass_G2(idx),STRADa_nuc_G2(idx),'poly1');
  h=plot(fit_line);
  h.Color = G2_color .*.7;
  h.LineWidth = 3;
  % Scatter plot by cell cycle MEAN
  h=scatter(cell_mass_EG1_mean,STRADa_nuc_EG1_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',EG1_color);
  h=scatter(cell_mass_LG1_mean,STRADa_nuc_LG1_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',LG1_color);
  h=scatter(cell_mass_S_mean,STRADa_nuc_S_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',S_color);
  h=scatter(cell_mass_G2_mean,STRADa_nuc_G2_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',G2_color);
  h=scatter(cell_mass_G2_mean,STRADa_nuc_G2_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',G2_color);

  legend({'EG1 Cell','LG1 Cell','S Cell','G2 Cell','EG1 Trend','LG1 Trend','S Trend','G2 Trend','EG1 Mean', 'LG1 Mean', 'S Mean', 'G2 Mean'});
  title('Cell Mass (SE) vs Nuclear STRADa (STRADa Normalized by SE)')
  xlabel('Cell Mass (SE)');
  ylabel('<- Less STRADa in Nucleus         Ratio         More STRADa in Nucleus ->');
  axis tight

  filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5 - binned plots/mass-vs-nuc.png');
  export_fig(filename, '-m2 -transparent -nocrop')

  
  
end