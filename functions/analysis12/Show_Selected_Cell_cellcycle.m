function fun(ResultTable, Data, resolution, special)
  channels = Data.O.Channels;
  STRADa_channel = find(ismember(channels, 'STRADa'));
  cell_mass_channel = find(ismember(channels, 'SE'));
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
  higher_prctile = 99;
  cell_size = set_outliers_to(cell_size, lower_prctile, higher_prctile, NaN);
  STRADa_nuc = set_outliers_to(STRADa_nuc, lower_prctile, higher_prctile, NaN);
  STRADa_total = set_outliers_to(STRADa_total, lower_prctile, higher_prctile, NaN);
  STRADa_cyto = set_outliers_to(STRADa_cyto, lower_prctile, higher_prctile, NaN);
  STRADa_localization_ratio_nuc_cyto = set_outliers_to(STRADa_localization_ratio_nuc_cyto, lower_prctile, higher_prctile, NaN);
  STRADa_localization_ratio_nuc_total = set_outliers_to(STRADa_localization_ratio_nuc_total, lower_prctile, higher_prctile, NaN);
  STRADa_localization_ratio_cyto_total = set_outliers_to(STRADa_localization_ratio_cyto_total, lower_prctile, higher_prctile, NaN);

  all_STRA = [STRADa_localization_ratio_nuc_cyto ...
              STRADa_localization_ratio_nuc_total ...
              STRADa_localization_ratio_cyto_total]; % used for consistent plot axis limits

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

  STRA_ratio_EG1_mean = nanmean(STRA_ratio_EG1);
  STRA_ratio_LG1_mean = nanmean(STRA_ratio_LG1);
  STRA_ratio_S_mean = nanmean(STRA_ratio_S);
  STRA_ratio_G2_mean = nanmean(STRA_ratio_G2);
  cell_size_EG1_mean = nanmean(cell_size_EG1);
  cell_size_LG1_mean = nanmean(cell_size_LG1);
  cell_size_S_mean = nanmean(cell_size_S);
  cell_size_G2_mean = nanmean(cell_size_G2);

  figure('Position',[251,    91,   883,   708]);
  hold on
  
  % Make a meshgrid for the heatmap
  % X = linspace(prctile(cell_size,1),prctile(cell_size,99),resolution);
  % Y = linspace(prctile(STRADa_localization_ratio_nuc_cyto,1),prctile(STRADa_localization_ratio_nuc_cyto,99),resolution);
  % [x, y] = meshgrid(X,Y);
  % % Plot Heat Map
  % two_dimensions = [ cell_size  STRADa_localization_ratio_nuc_cyto ];
  % [d, f]= ksdensity(two_dimensions, [x(:) y(:)]);
  % JPD = reshape(d,[resolution resolution]); % JDP stands for joint_probability_density
  % if strcmp(special, 'log')
  %   JPD = log(JPD);
  % end
%  h=imagesc(X,Y,JPD);colormap(gca,jet)
  %alpha(h,.5);

  % Scatter plot by cell cycle
  X=cell_size_EG1(1:50:end);
  Y=STRA_ratio_EG1(1:50:end);
  h=scatter(X,Y,5,'filled', 'MarkerFaceAlpha',1,'MarkerFaceColor',[0 .6 0]);
  X=cell_size_LG1(1:50:end);
  Y=STRA_ratio_LG1(1:50:end);
  h=scatter(X,Y,5,'filled', 'MarkerFaceAlpha',1,'MarkerFaceColor',[0 1 .5]);
  X=cell_size_S(1:50:end);
  Y=STRA_ratio_S(1:50:end);
  h=scatter(X,Y,5,'filled', 'MarkerFaceAlpha',1,'MarkerFaceColor',[.2 1 1]);
  X=cell_size_G2(1:50:end);
  Y=STRA_ratio_G2(1:50:end);
  h=scatter(X,Y,5,'filled', 'MarkerFaceAlpha',1,'MarkerFaceColor',[.4 .6 1]);
  % Scatter plot by cell cycle MEAN
  h=scatter(cell_size_EG1_mean,STRA_ratio_EG1_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',[0 .6 0]);
  h=scatter(cell_size_LG1_mean,STRA_ratio_LG1_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',[0 1 .5]);
  h=scatter(cell_size_S_mean,STRA_ratio_S_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',[.2 1 1]);
  h=scatter(cell_size_G2_mean,STRA_ratio_G2_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',[.4 .6 1]);
  % Make plot pretty
  legend({'EG1', 'LG1', 'S', 'G2','EG1 Mean', 'LG1 Mean', 'S Mean', 'G2 Mean'});
  set(gca,'ydir','normal')
  title('STRADa Localization Ratio (average intensity) vs Cell Area')
  xlabel('Cell Size (area)');ylabel('STRADa_n_u_c/STRADa_c_y_t_o');
  legendmarkeradjust(90);


  % Choose cell of interest by clicking data point in plot
  [x,y,button] = ginput(1);

  X = cell_size;
  Y = STRADa_localization_ratio_nuc_cyto;

  % Calc closest data point to click
  closest = Inf;
  closest_id = NaN;
  for i=1:height(ResultTable)
    x_difference = abs(x - X(i));
    y_difference = abs(y - Y(i));
    total_difference = x_difference + y_difference;
    if total_difference < closest
      closest = total_difference;
      closest_id = i;
    end
  end

  closest_x = ResultTable.Centroid(closest_id,1);
  closest_y = ResultTable.Centroid(closest_id,2);


  %% Find and segment single image.
  ImageIDs=Data.O.ImageIDs;
  t1=datetime;
  t2=datetime;

  % Get image this cell is in
  I= ImageIDs.Row==ResultTable.Row(closest_id) ...
      & ImageIDs.Column==ResultTable.Column(closest_id) ...
      & ImageIDs.Field==ResultTable.Field(closest_id) ...
      & ImageIDs.Channel==1;

  addpath('Z:\OPRETTA\Operetta Image Processing\User Interface\lib\NewVersion')
  [iterTable,ImageID,O]=O_SegmentCells_v6_SingleImage(Data.O,ImageIDs(I,:),t1,t2);
  rmpath('Z:\OPRETTA\Operetta Image Processing\User Interface\lib\NewVersion')
  
  ImageID.FileName

  %% Display each channel with an x to mark the cell and nuc and cyto boundries
  set(0,'DefaultFigureWindowStyle','docked');
  warning off images:initSize:adjustingMag;
  warning off images:imshow:magnificationMustBeFitForDockedFigure;
  nuc_perim = bwperim(O.BW{DAPI_channel});
  cyto_perim = bwperim(O.BW{cell_mass_channel});

  figure(137)
  for cid=1:length(channels)
    subplot(2,2,cid)
    % display this channel
    img=O.IM{cid};
    imshow(img, [prctile(img(:), 3), prctile(img(:), 99.9)]);
    title(channels(cid))
    hold on
    % nuc boundries
    [Y,X] = find(nuc_perim);
    h=scatter(X,Y,4,'filled','MarkerFaceAlpha',3/4);
    % cyto boundries
    [Y,X] = find(cyto_perim);
    h=scatter(X,Y,4,'filled','MarkerFaceAlpha',2/3);
    % mark the x

    text(closest_x-5,closest_y-5,'x','Color','r')

    xlim([closest_x - 120 closest_x + 120])
    ylim([closest_y - 120 closest_y + 120])

  end
end

