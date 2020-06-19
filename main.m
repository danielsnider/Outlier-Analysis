addpath(genpath('functions'));

results_with_PxId_and_prcile = '\\carbon.research.sickkids.ca\rkafri\OPRETTA\Operetta Processed OutPutFiles\Dataset_20170222_STRADa_SE_Eden_with_prctile_and_PxIdRESULTS\ResultTable2.mat';
%load('ResultTable.mat')
load(results_with_PxId_and_prcile)
load('SegPar_SegmentationParameters.mat')


%% Analysis 1 and 2
% Basics
basic_stats(ResultTable,Data);
display_sample_images(ResultTable,Data);

% Single variable analysis
each_channel_ksdensity(ResultTable,Data); % for nuc and cyto
STRADa_ksdensity(ResultTable,Data); % and cell size ksdensity

% Cell mass vs STRADa
cell_mass_vs_STRADa_localization_difference_scatter(ResultTable,Data);
cell_mass_vs_STRADa_scaled_for_size_contour(ResultTable,Data,20);
cell_mass_vs_STRADa_scaled_for_size_JPD(ResultTable,Data,20);
cell_cycle_JPD_nuc_cyto_cell_mass(ResultTable,Data,15);

% Cell size
cellsize_vs_STRADa_localization_ratio_nuc_cyto_JPD(ResultTable,Data,90,''); %110
cellsize_vs_STRADa_localization_ratio_3_subplots_JPD(ResultTable,Data,20,''); % 90
cell_cycle_nuc_cyto_cellsize_JPD(ResultTable,Data,10,''); % 60

cellsize_vs_STRADa_localization_ratio_nuc_cyto_JPD_AND_scatter(ResultTable,Data,20,'');


cellsize_vs_STRADa_scaled_for_size_contour(ResultTable,Data,75);
cellsize_vs_STRADa_nuc_cyto_scatter(ResultTable,Data);
cellsize_vs_STRADa_scaled_nuc_cyto_scatter(ResultTable,Data);
cellsize_vs_STRADa_nuc_vs_STRADa_cyto_scatter3(ResultTable,Data);

% Find single cell images using a plot (interactive)
Show_Selected_Cell(ResultTable,Data);
Show_Selected_Cell_cellcycle(ResultTable,Data,10,'');

%% Analysis 3
%single_cell_overview(ResultTable,Data)
single_cell_overview_by_mass(ResultTable,Data)

%% Analysis 4
small_large_vs_nuc_cyto_STRADa(ResultTable, Data, 10, '')

%% Analysis 5
mass_vs_STRADa(ResultTable, Data)
mass_vs_STRADa_binned(ResultTable, Data)

%% Analysis 6
mass_vs_STRADa_a6(ResultTable, Data)
mass_vs_STRADa_binned_a6(ResultTable, Data)
single_cell_overview_by_ratio(ResultTable,Data)
single_cell_overview_top_and_bottom(ResultTable,Data)

%% Analysis 7 - prctile
single_cell_overview_prctile(ResultTable,Data)