function fun(ResultTable, Data, resolution)
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


  % Get stats for each cell stage
  cell_size_EG1 = cell_size(idx_EG1);
  cell_size_LG1 = cell_size(idx_LG1);
  cell_size_S = cell_size(idx_S);
  cell_size_G2 = cell_size(idx_G2);
  STRA_total_EG1 = STRADa_total(idx_EG1);
  STRA_total_LG1 = STRADa_total(idx_LG1);
  STRA_total_S = STRADa_total(idx_S);
  STRA_total_G2 = STRADa_total(idx_G2);
  STRA_nuc_EG1 = STRADa_nuc(idx_EG1);
  STRA_nuc_LG1 = STRADa_nuc(idx_LG1);
  STRA_nuc_S = STRADa_nuc(idx_S);
  STRA_nuc_G2 = STRADa_nuc(idx_G2);
  STRA_cyto_EG1 = STRADa_cyto(idx_EG1);
  STRA_cyto_LG1 = STRADa_cyto(idx_LG1);
  STRA_cyto_S = STRADa_cyto(idx_S);
  STRA_cyto_G2 = STRADa_cyto(idx_G2);
  STRA_ratio_EG1 = STRADa_localization_ratio_nuc_cyto(idx_EG1);
  STRA_ratio_LG1 = STRADa_localization_ratio_nuc_cyto(idx_LG1);
  STRA_ratio_S = STRADa_localization_ratio_nuc_cyto(idx_S);
  STRA_ratio_G2 = STRADa_localization_ratio_nuc_cyto(idx_G2);
  

  figure('Position',[251,    91,   1083,   708]);
  hold on
  ksdensity(cell_size_EG1)
  ksdensity(cell_size_LG1)
  ksdensity(cell_size_S)
  ksdensity(cell_size_G2)
  title('Cell Area ksdensity')
  xlabel('Cell Area')
  ylabel('Probability')
  legend({'Cell Area in EG1', 'Cell Area in LG1', 'Cell Area in S', 'Cell Area in G2'})
  
  % Save the plot to disk
  filename = sprintf('cell_size_ksdensity.png');
  export_fig(['plots\' filename], '-m2')

  figure('Position',[251,    91,   1083,   708]);
  subplot(2,2,1)
  hold on
  ksdensity(STRA_total_EG1)
  ksdensity(STRA_total_LG1)
  ksdensity(STRA_total_S)
  ksdensity(STRA_total_G2)
  title('Total STRADa ksdensity')
  xlabel('Total STRADa avg int')
  ylabel('Probability')
  legend({'Total STRADa avg int (EG1)', 'Total STRADa avg int (LG1)', 'Total STRADa avg int (S)', 'Total STRADa avg int (G2)'})

  subplot(2,2,2)
  hold on
  ksdensity(STRA_ratio_EG1)
  ksdensity(STRA_ratio_LG1)
  ksdensity(STRA_ratio_S)
  ksdensity(STRA_ratio_G2)
  title('Ratio (nuc/cyto) ksdensity')
  xlabel('Ratio (nuc/cyto) STRADa')
  ylabel('Probability')
  legend({'Ratio (nuc/cyto) STRADa avg int (EG1)', 'Ratio (nuc/cyto) STRADa avg int (LG1)', 'Ratio (nuc/cyto) STRADa avg int (S)', 'Ratio (nuc/cyto) STRADa avg int (G2)'})

  subplot(2,2,3)
  hold on
  ksdensity(STRA_nuc_EG1)
  ksdensity(STRA_nuc_LG1)
  ksdensity(STRA_nuc_S)
  ksdensity(STRA_nuc_G2)
  title('Nuclear STRADa ksdensity')
  xlabel('Nuclear STRADa avg int')
  ylabel('Probability')
  legend({'Nuclear STRADa avg int (EG1)', 'Nuclear STRADa avg int (LG1)', 'Nuclear STRADa avg int (S)', 'Nuclear STRADa avg int (G2)'})

  subplot(2,2,4)
  hold on
  ksdensity(STRA_cyto_EG1)
  ksdensity(STRA_cyto_LG1)
  ksdensity(STRA_cyto_S)
  ksdensity(STRA_cyto_G2)
  title('Cytoplasmic STRADa ksdensity')
  xlabel('Cytoplasmic STRADa avg int')
  ylabel('Probability')
  legend({'Cytoplasmic STRADa avg int (EG1)', 'Cytoplasmic STRADa avg int (LG1)', 'Cytoplasmic STRADa avg int (S)', 'Cytoplasmic STRADa avg int (G2)'})


  % Save the plot to disk
  filename = sprintf('cell_cycle_ksdensity.png');
  export_fig(['plots\' filename], '-m2')
end
