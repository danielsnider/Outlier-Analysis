 function fun(ResultTable, Data)
  channels = Data.O.Channels;
  cell_stages = {'EG1', 'LG1', 'S', 'G2'};
  STRADa_channel = find(ismember(channels, 'STRADa'));
  cell_mass_channel = find(ismember(channels, 'SE'));
  DAPI_channel = find(ismember(channels, 'DAPI'));
  geminin_channel = find(ismember(channels, 'Geminin'));
  EG1_color = [0/255 255/255 150/255];
  LG1_color = [232/255 227/255 12/255];
  S_color = [255/255 102/255 0/255];
  G2_color = [207/255 12/255 232/255];
  DAPI = ResultTable.NInt(:,DAPI_channel);
  geminin = ResultTable.NInt(:,geminin_channel);
  cell_mass = ResultTable.CInt(:,cell_mass_channel);
  nuc_mass = ResultTable.NInt(:,cell_mass_channel);
  cyto_mass = cell_mass - nuc_mass;
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;
  STRADa_nuc = STRADa_nuc./nuc_mass;
  STRADa_total = STRADa_total./cell_mass;
  STRADa_cyto = STRADa_cyto./cyto_mass;
  STRADa_ratio = STRADa_nuc./STRADa_cyto;
  STRADa_ratio_nuc_total = STRADa_nuc./STRADa_total;
  STRADa_ratio_cyto_total = STRADa_cyto./STRADa_total;

  % % Eliminate outliers
  lower_prctile = 0.1;
  higher_prctile = 99.9;
  cell_mass = set_outliers_to(cell_mass, lower_prctile, higher_prctile, NaN);
  nuc_mass = set_outliers_to(nuc_mass, lower_prctile, higher_prctile, NaN);
  cyto_mass = set_outliers_to(cyto_mass, lower_prctile, higher_prctile, NaN);
  STRADa_nuc = set_outliers_to(STRADa_nuc, lower_prctile, higher_prctile, NaN);
  STRADa_total = set_outliers_to(STRADa_total, lower_prctile, higher_prctile, NaN);
  STRADa_cyto = set_outliers_to(STRADa_cyto, lower_prctile, higher_prctile, NaN);
  STRADa_nuc = set_outliers_to(STRADa_nuc, lower_prctile, higher_prctile, NaN);
  STRADa_total = set_outliers_to(STRADa_total, lower_prctile, higher_prctile, NaN);
  STRADa_cyto = set_outliers_to(STRADa_cyto, lower_prctile, higher_prctile, NaN);
  STRADa_ratio = set_outliers_to(STRADa_ratio, lower_prctile, higher_prctile, NaN);
  STRADa_ratio_nuc_total = set_outliers_to(STRADa_ratio_nuc_total, lower_prctile, higher_prctile, NaN);
  STRADa_ratio_cyto_total = set_outliers_to(STRADa_ratio_cyto_total, lower_prctile, higher_prctile, NaN);

  

  % Find cell stages
  [idx_EG1,idx_LG1,idx_G1S,idx_S,idx_G2] = FindStages_v2(DAPI,log(geminin));

  % combine LG1 and G1S into one group (LG1) because G1S wasn't detected well enough
  idx_LG1 = idx_LG1 | idx_G1S;
  idx_G1S = []; % delete so no mistake is made


  % Get stats for each cell stage
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
  STRADa_ratio_EG1 = STRADa_ratio(idx_EG1);
  STRADa_ratio_LG1 = STRADa_ratio(idx_LG1);
  STRADa_ratio_S = STRADa_ratio(idx_S);
  STRADa_ratio_G2 = STRADa_ratio(idx_G2);


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






  n=1

  %List of aspects to be visualized
  aspects(n).num_bins = 7;
  aspects(n).X = 'STRADa_ratio';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = true;
  aspects(n).title = 'Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).xlabel = '<- More STRADa in Cytoplasm                          Ratio                          More STRADa in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 29;
  aspects(n).X = 'STRADa_ratio';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = true;
  aspects(n).title = 'Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).xlabel = '<- More STRADa in Cytoplasm                          Ratio                          More STRADa in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 129;
  aspects(n).X = 'STRADa_ratio';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = true;
  aspects(n).title = 'Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).xlabel = '<- More STRADa in Cytoplasm                          Ratio                          More STRADa in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 7;
  aspects(n).X = 'STRADa_ratio';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = false;
  aspects(n).title = 'Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).xlabel = '<- More STRADa in Cytoplasm                          Ratio                          More STRADa in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 29;
  aspects(n).X = 'STRADa_ratio';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = false;
  aspects(n).title = 'Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).xlabel = '<- More STRADa in Cytoplasm                          Ratio                          More STRADa in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 129;
  aspects(n).X = 'STRADa_ratio';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = false;
  aspects(n).title = 'Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).xlabel = '<- More STRADa in Cytoplasm                          Ratio                          More STRADa in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 7;
  aspects(n).X = 'cell_mass';
  aspects(n).Y = 'STRADa_ratio';
  aspects(n).hide_points = true;
  aspects(n).title = 'Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)';
  aspects(n).xlabel = 'Cell Mass (SE)';
  aspects(n).ylabel = '<- More STRADa in Cytoplasm                          Ratio                          More STRADa in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 29;
  aspects(n).X = 'cell_mass';
  aspects(n).Y = 'STRADa_ratio';
  aspects(n).hide_points = true;
  aspects(n).title = 'Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)';
  aspects(n).xlabel = 'Cell Mass (SE)';
  aspects(n).ylabel = '<- More STRADa in Cytoplasm                          Ratio                          More STRADa in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 129;
  aspects(n).X = 'cell_mass';
  aspects(n).Y = 'STRADa_ratio';
  aspects(n).hide_points = true;
  aspects(n).title = 'Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)';
  aspects(n).xlabel = 'Cell Mass (SE)';
  aspects(n).ylabel = '<- More STRADa in Cytoplasm                          Ratio                          More STRADa in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 7;
  aspects(n).X = 'cell_mass';
  aspects(n).Y = 'STRADa_ratio';
  aspects(n).hide_points = false;
  aspects(n).title = 'Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)';
  aspects(n).xlabel = 'Cell Mass (SE)';
  aspects(n).ylabel = '<- More STRADa in Cytoplasm                          Ratio                          More STRADa in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 29;
  aspects(n).X = 'cell_mass';
  aspects(n).Y = 'STRADa_ratio';
  aspects(n).hide_points = false;
  aspects(n).title = 'Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)';
  aspects(n).xlabel = 'Cell Mass (SE)';
  aspects(n).ylabel = '<- More STRADa in Cytoplasm                          Ratio                          More STRADa in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 129;
  aspects(n).X = 'cell_mass';
  aspects(n).Y = 'STRADa_ratio';
  aspects(n).hide_points = false;
  aspects(n).title = 'Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)';
  aspects(n).xlabel = 'Cell Mass (SE)';
  aspects(n).ylabel = '<- More STRADa in Cytoplasm                          Ratio                          More STRADa in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 7;
  aspects(n).X = 'STRADa_nuc';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = true;
  aspects(n).title = 'Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).xlabel = '<- Less in Nucleus          Nuclear STRADa          More in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 29;
  aspects(n).X = 'STRADa_nuc';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = true;
  aspects(n).title = 'Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).xlabel = '<- Less in Nucleus          Nuclear STRADa          More in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 129;
  aspects(n).X = 'STRADa_nuc';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = true;
  aspects(n).title = 'Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).xlabel = '<- Less in Nucleus          Nuclear STRADa          More in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 7;
  aspects(n).X = 'STRADa_nuc';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = false;
  aspects(n).title = 'Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).xlabel = '<- Less in Nucleus          Nuclear STRADa          More in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 29;
  aspects(n).X = 'STRADa_nuc';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = false;
  aspects(n).title = 'Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).xlabel = '<- Less in Nucleus          Nuclear STRADa          More in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 129;
  aspects(n).X = 'STRADa_nuc';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = false;
  aspects(n).title = 'Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).xlabel = '<- Less in Nucleus          Nuclear STRADa          More in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 7;
  aspects(n).X = 'cell_mass';
  aspects(n).Y = 'STRADa_nuc';
  aspects(n).hide_points = true;
  aspects(n).title = 'Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)';
  aspects(n).xlabel = '<- Less in Nucleus          Nuclear STRADa          More in Nucleus ->';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 29;
  aspects(n).X = 'cell_mass';
  aspects(n).Y = 'STRADa_nuc';
  aspects(n).hide_points = true;
  aspects(n).title = 'Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)';
  aspects(n).xlabel = '<- Less in Nucleus          Nuclear STRADa          More in Nucleus ->';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 129;
  aspects(n).X = 'cell_mass';
  aspects(n).Y = 'STRADa_nuc';
  aspects(n).hide_points = true;
  aspects(n).title = 'Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)';
  aspects(n).xlabel = '<- Less in Nucleus          Nuclear STRADa          More in Nucleus ->';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 7;
  aspects(n).X = 'cell_mass';
  aspects(n).Y = 'STRADa_nuc';
  aspects(n).hide_points = false;
  aspects(n).title = 'Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)';
  aspects(n).xlabel = '<- Less in Nucleus          Nuclear STRADa          More in Nucleus ->';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 29;
  aspects(n).X = 'cell_mass';
  aspects(n).Y = 'STRADa_nuc';
  aspects(n).hide_points = false;
  aspects(n).title = 'Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)';
  aspects(n).xlabel = '<- Less in Nucleus          Nuclear STRADa          More in Nucleus ->';
  aspects(n).ylabel = 'Cell Mass (SE)';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;

  aspects(n).num_bins = 129;
  aspects(n).X = 'STRADa_nuc';
  aspects(n).Y = 'cell_mass';
  aspects(n).hide_points = false;
  aspects(n).title = 'Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)';
  aspects(n).xlabel = 'Cell Mass (SE)';
  aspects(n).ylabel = '<- Less in Nucleus          Nuclear STRADa          More in Nucleus ->';
  aspects(n).filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/%s_vs_%s_bins%s_hidepoints_%s.png', aspects(n).X,aspects(n).Y,num2str(aspects(n).num_bins),string(aspects(n).hide_points));
  aspects(n).resolution = [2          42        1191        1074];
  n=n+1;
































  for i=1:size(aspects,2)
    num_bins = aspects(i).num_bins;
    num_cell_stages = 4;
    EG1_means_per_bin = nan(num_bins,1);
    LG1_means_per_bin = nan(num_bins,1);
    S_means_per_bin = nan(num_bins,1);
    G2_means_per_bin = nan(num_bins,1);
    EG1_errors_per_bin = nan(num_bins,1);
    LG1_errors_per_bin = nan(num_bins,1);
    S_errors_per_bin = nan(num_bins,1);
    G2_errors_per_bin = nan(num_bins,1);
    bin_centers = nan(num_bins,1);
    % Error = std / sqrt(length(data_points))
    num_columns = num_bins*num_cell_stages;
    data_for_bars = nan(length(cell_mass),num_columns); % data for bar chart is represented as columns for each bar and rows for each datum. Since some bars will have more data than others, some columns will have a lot of NaNs (to keep the matrix shape rectangular)
    for bin_num=1:num_bins
        bin_num
      % Calc bin spacing
      X = eval(aspects(i).X);
      Y = eval(aspects(i).Y);
      bin_spacing = (max(X) - min(X))/num_bins;
      bin_start = min(X) + (bin_spacing * (bin_num-1))
      bin_end = min(X) + (bin_spacing * bin_num);
      bin_centers(bin_num) = bin_end - bin_spacing/2;
      % Find cells in bin by cell stage
      in_bin = bin_start < X & X < bin_end;
      EG1_in_bin = idx_EG1 & in_bin;
      LG1_in_bin = idx_LG1 & in_bin;
      S_in_bin = idx_S & in_bin;
      G2_in_bin = idx_G2 & in_bin;
      % Get cell mass
      value_EG1_in_bin = Y(EG1_in_bin);
      value_LG1_in_bin = Y(LG1_in_bin);
      value_S_in_bin = Y(S_in_bin);
      value_G2_in_bin = Y(G2_in_bin);
      % Store data points appropriately
      % Four columns will be filled with whatever amount of data is present for the cell stage and this bin
      data_for_bars(1:length(value_EG1_in_bin),((bin_num-1)*num_cell_stages)+1) = value_EG1_in_bin;
      data_for_bars(1:length(value_LG1_in_bin),((bin_num-1)*num_cell_stages)+2) = value_LG1_in_bin;
      data_for_bars(1:length(value_S_in_bin),((bin_num-1)*num_cell_stages)+3) = value_S_in_bin;
      data_for_bars(1:length(value_G2_in_bin),((bin_num-1)*num_cell_stages)+4) = value_G2_in_bin;
      % Store means
      EG1_means_per_bin(bin_num) = nanmean(value_EG1_in_bin);
      LG1_means_per_bin(bin_num) = nanmean(value_LG1_in_bin);
      S_means_per_bin(bin_num) = nanmean(value_S_in_bin);
      G2_means_per_bin(bin_num) = nanmean(value_G2_in_bin);
      % Store errors
      EG1_errors_per_bin(bin_num) = nanstd(value_EG1_in_bin)/sqrt(length(value_EG1_in_bin));
      LG1_errors_per_bin(bin_num) = nanstd(value_LG1_in_bin)/sqrt(length(value_LG1_in_bin));
      S_errors_per_bin(bin_num) = nanstd(value_S_in_bin)/sqrt(length(value_S_in_bin));
      G2_errors_per_bin(bin_num) = nanstd(value_G2_in_bin)/sqrt(length(value_G2_in_bin));
      %%% Plot to debug bins
      % x = [bin_start bin_start];
      % y = [0 1];
      % line(x,y,'Color','blue');
      % % pause
      % x = [bin_end bin_end];
      % y = [0 1];
      % line(x,y,'Color','red');
      % % pause
    end

    num_plot_columns = num_bins*(num_cell_stages+1);
    spacing = 1:num_plot_columns;
    for num=4:5:num_columns*(num_cell_stages+1)
       spacing(spacing==num+1) = [];
    end
    figure('Position',aspects(i).resolution);
    % figure('Position',[1          41        1920        1083]);
    H=notBoxPlot(data_for_bars,spacing,'jitter',0.6);
    d=[H.data];
    hold on

    % Disable the displaying the data points
    if aspects(i).hide_points
      set(d,'Visible','off')
    end


    % Box Colors
    EG1_IND=zeros(1,num_columns);
    EG1_IND(1:4:num_columns)=1;
    LG1_IND=zeros(1,num_columns);
    LG1_IND(2:4:num_columns)=1;
    S_IND=zeros(1,num_columns);
    S_IND(3:4:num_columns)=1;
    G2_IND=zeros(1,num_columns);
    G2_IND(4:4:num_columns)=1;

    set([H(find(EG1_IND)).data],'MarkerSize',4,...
        'markerFaceColor',EG1_color*0.5,...
        'markerEdgeColor', 'none')
    set([H(find(EG1_IND)).semPtch],...
        'FaceColor',EG1_color*0.25,...
        'EdgeColor','none')
    set([H(find(EG1_IND)).sdPtch],...
        'FaceColor',EG1_color*0.75,...
        'EdgeColor','none')
    set([H(find(EG1_IND)).mu],...
        'Color','b')

    set([H(find(LG1_IND)).data],'MarkerSize',4,...
        'markerFaceColor',LG1_color*0.5,...
        'markerEdgeColor', 'none')
    set([H(find(LG1_IND)).semPtch],...
        'FaceColor',LG1_color*0.25,...
        'EdgeColor','none')
    set([H(find(LG1_IND)).sdPtch],...
        'FaceColor',LG1_color*0.75,...
        'EdgeColor','none')
    set([H(find(LG1_IND)).mu],...
        'Color','b')

    set([H(find(S_IND)).data],'MarkerSize',4,...
        'markerFaceColor',S_color*0.5,...
        'markerEdgeColor', 'none')
    set([H(find(S_IND)).semPtch],...
        'FaceColor',S_color*0.25,...
        'EdgeColor','none')
    set([H(find(S_IND)).sdPtch],...
        'FaceColor',S_color*0.75,...
        'EdgeColor','none')
    set([H(find(S_IND)).mu],...
        'Color','b')

    set([H(find(G2_IND)).data],'MarkerSize',4,...
        'markerFaceColor',G2_color*0.5,...
        'markerEdgeColor', 'none')
    set([H(find(G2_IND)).semPtch],...
        'FaceColor',G2_color*0.25,...
        'EdgeColor','none')
    set([H(find(G2_IND)).sdPtch],...
        'FaceColor',G2_color*0.75,...
        'EdgeColor','none')
    set([H(find(G2_IND)).mu],...
        'Color','b')

    % Trend lines
    bin_centers_on_plot  = (2:num_cell_stages+1:num_plot_columns)+.5;
    plot(bin_centers_on_plot,EG1_means_per_bin,'color',EG1_color*.9,'LineWidth', 2);
    plot(bin_centers_on_plot,LG1_means_per_bin,'color',LG1_color*.9,'LineWidth', 2);
    plot(bin_centers_on_plot,S_means_per_bin,'color',S_color*.9,'LineWidth', 2);
    plot(bin_centers_on_plot,G2_means_per_bin,'color',G2_color*.9,'LineWidth', 2);

    % Mean line color
    set([H.mu],'color',[0.5 .5 .5])

    % X ticks
    xticklabels = {};
    for bin_num=1:round(num_bins/8):num_bins
      if max(bin_centers) > 10
        xticklabel = sprintf('%.1e', bin_centers(bin_num));
      else
        xticklabel = sprintf('%.4f', bin_centers(bin_num));
      end
      % xticklabel = sprintf('%d',bin_num);
      xticklabels{bin_num} = xticklabel;
    end
    set(gca,'XTick',bin_centers_on_plot) % The middle point of each bin
    set(gca,'XTickLabels',xticklabels)


    % Pretty
    box on
    title(aspects(i).title);
    ylabel(aspects(i).ylabel);
    xlabel(aspects(i).xlabel);

    pause(.2) % wait to load figure
    export_fig(aspects(i).filename, '-m2 -transparent')
%    close all
  end


  
end

