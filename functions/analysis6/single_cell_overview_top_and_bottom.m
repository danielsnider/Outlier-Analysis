function fun(ResultTable, Data)
  %% Filter Errors
  % Remove if cyto area and nuc area is the same
  ResultTable(ResultTable.CArea == ResultTable.NArea,:) = [];

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
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;



  % SET NUMBER OF CELLS TO LOOK AT
  num_cells = length(STRADa_cyto); % how many cells do you want pictures for? double this number because we'll look at the upper and lower outliers of this many cells
  multiplier = 1; % how many jumps between sequential cells to look at
  % Get n number of cells of the upper and lower outliers

  n=1
  %List of aspects to be visualized
  aspects(n).code = 'STRADa_ratio = STRADa_nuc./STRADa_total;';
  aspects(n).foldername = 'STRADa_ratio = STRADa_nuc DB STRADa_total';
  n=n+1;

  aspects(n).code = 'STRADa_ratio = (STRADa_nuc./nuc_size) ./ (STRADa_cyto./cyto_size);';
  aspects(n).foldername = 'STRADa_ratio = (STRADa_nuc DB nuc_size) DB (STRADa_cyto DB cyto_size)';
  n=n+1;

  aspects(n).code = 'STRADa_ratio = STRADa_nuc;';
  aspects(n).foldername = 'STRADa_ratio = STRADa_nuc';
  n=n+1;

  aspects(n).code = 'STRADa_ratio = STRADa_nuc./nuc_size;';
  aspects(n).foldername = 'STRADa_ratio = STRADa_nuc DB nuc_size';
  n=n+1;

  aspects(n).code = 'STRADa_ratio = (STRADa_nuc./nuc_mass) ./ (STRADa_cyto./cyto_mass);';
  aspects(n).foldername = 'STRADa_ratio = (STRADa_nuc DB nuc_mass) DB (STRADa_cyto DB cyto_mass)';
  n=n+1;

  figure('Position',[251,    91,   883,   708]);
  for i=1:size(aspects,2)
    eval(aspects(i).code)

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
    STRADa_ratio_EG1 = STRADa_ratio(idx_EG1);
    STRADa_ratio_LG1 = STRADa_ratio(idx_LG1);
    STRADa_ratio_S = STRADa_ratio(idx_S);
    STRADa_ratio_G2 = STRADa_ratio(idx_G2);


    % % Eliminate outliers
    lower_prctile = 1;
    higher_prctile = 99;
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
    STRADa_ratio = set_outliers_to(STRADa_ratio, lower_prctile, higher_prctile, NaN);


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


    % SORT TABLE BY RATIO
    [sorted_data,index] = sortrows(STRADa_ratio);
    ids_of_lowest = index(1:multiplier:num_cells*multiplier);
    ids_of_highest = index(length(index)-(num_cells*multiplier)+1:multiplier:length(index));
    cell_ids = [ids_of_highest; ids_of_lowest];
    max_count = 10;
    low_count = 0;
    high_count = 0;

    for c=1:length(cell_ids)
      if mod(c,2) == 0 %number is even
        cell_id = cell_ids(c);
      else
        cell_id = cell_ids(end-c);
      end
      cell_stage = cell_stages(find([idx_EG1(cell_id) idx_LG1(cell_id) idx_S(cell_id) idx_G2(cell_id)]));
      if isempty(cell_stage)
          continue
      end
      if isnan(STRADa_ratio(cell_id))
        continue
      end
      if isnan(cell_mass(cell_id))
        continue
      end
      if mod(c,2) == 0 %number is even
        if high_count >= max_count
          'skip high cell'
          continue
        end
      else
        if low_count >= max_count
          'skip low cell'
          continue
        end
      end

      clf
      subplot(3,3,5)
      hold on
      % Scatter plot by cell cycle
      idx = isfinite(cell_mass_EG1) & isfinite(STRADa_ratio_EG1);
      X=cell_mass_EG1(idx);
      Y=STRADa_ratio_EG1(idx);
      h=plot(X(1:50:end),Y(1:50:end),'.','color',EG1_color);
      idx = isfinite(cell_mass_LG1) & isfinite(STRADa_ratio_LG1);
      X=cell_mass_LG1(idx);
      Y=STRADa_ratio_LG1(idx);
      h=plot(X(1:50:end),Y(1:50:end),'.','color',LG1_color);
      idx = isfinite(cell_mass_S) & isfinite(STRADa_ratio_S);
      X=cell_mass_S(idx);
      Y=STRADa_ratio_S(idx);
      h=plot(X(1:50:end),Y(1:50:end),'.','color',S_color);
      idx = isfinite(cell_mass_G2) & isfinite(STRADa_ratio_G2);
      X=cell_mass_G2(idx);
      Y=STRADa_ratio_G2(idx);
      h=plot(X(1:50:end),Y(1:50:end),'.','color',G2_color);
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
      title('Cell Mass vs STRADa localization')
      xlabel('Cell Mass (SE)');
      ylabel('% of STRADa in Nucleus');
      xlim([min(cell_mass) max(cell_mass)])
      ylim([min(STRADa_ratio) max(STRADa_ratio)])
      
      % Plot an X for this cell
      X = cell_mass(cell_id);
      Y = STRADa_ratio(cell_id);
      text(X,Y,'x','Color','r')
      text(X,Y,'o','Color','k')
      grid on

      centroid_x = ResultTable.Centroid(cell_id,1);
      centroid_y = ResultTable.Centroid(cell_id,2);

      %% Find and segment single image.
      ImageIDs=Data.O.ImageIDs;
      t1=datetime;
      t2=datetime;

      % Get image this cell is in
      I= ImageIDs.Row==ResultTable.Row(cell_id) ...
          & ImageIDs.Column==ResultTable.Column(cell_id) ...
          & ImageIDs.Field==ResultTable.Field(cell_id) ...
          & ImageIDs.Channel==1;

      addpath('Z:\OPRETTA\Operetta Image Processing\User Interface\lib\NewVersion')
      [iterTable,ImageID,O]=O_SegmentCells_v6_SingleImage(Data.O,ImageIDs(I,:),t1,t2);
      rmpath('Z:\OPRETTA\Operetta Image Processing\User Interface\lib\NewVersion')
      
      nuc_perim = bwperim(O.BW{DAPI_channel});
      cyto_perim = bwperim(O.BW{cell_mass_channel});
      
      cyto_mask = O.BW{cell_mass_channel};
      labelled_id = cyto_mask(round(centroid_y), round(centroid_x));
      cyto_mask = cyto_mask==labelled_id;

      %% Display each channel with an x to mark the cell and nuc and cyto boundries
      set(0,'DefaultFigureWindowStyle','docked');
      warning off images:initSize:adjustingMag;
      warning off images:imshow:magnificationMustBeFitForDockedFigure;

      for cid=1:length(channels)
        subplot(3,3,cid)
        % display this channel
        img=O.IM{cid};
        if cid == cell_mass_channel
          imshow(img, [prctile(img(:), 3), prctile(img(:), 98)]);
        else
          imshow(img, [prctile(img(:), 3), prctile(img(:), 99.9)]);
        end
        title(channels(cid))
        hold on
        % nuc boundries
        [Y,X] = find(nuc_perim);
        h=scatter(X,Y,4,'filled','MarkerFaceAlpha',3/4);
        % cyto boundries
        [Y,X] = find(cyto_perim);
        h=scatter(X,Y,4,'filled','MarkerFaceAlpha',2/3);
        % mark the x

        text(centroid_x-1.00,centroid_y-1.00,'x','Color','r')

        xlim([centroid_x - 120 centroid_x + 120])
        ylim([centroid_y - 120 centroid_y + 120])

      end

      % Image of STRADa normalized by SE
      subplot(3,3,7)
      STRADa_img=O.IM{STRADa_channel};
      STRADa_img(find(~cyto_mask))=0;
      cell_mass_img=O.IM{cell_mass_channel};
      cell_mass_img(~cyto_mask)=0;
      img = STRADa_img ./ cell_mass_img;
      imshow(img, [prctile(img(:), 3), prctile(img(:), 97)]);
      title('STRADa Normalized by SE (STRADa./SE)')
      hold on
      % nuc boundries
      [Y,X] = find(nuc_perim);
      h=scatter(X,Y,4,'filled','MarkerFaceAlpha',3/4);
      % cyto boundries
      [Y,X] = find(cyto_perim);
      h=scatter(X,Y,4,'filled','MarkerFaceAlpha',2/3);
      xlim([centroid_x - 120 centroid_x + 120])
      ylim([centroid_y - 120 centroid_y + 120])

      %% Extra details
      subplot(3,3,8)
      hold on
      text(0,.5,['Cell Stage: ' char(cell_stage)]); 
      text(0,.4,['Cell Size Percentile (Overall): ' num2str(invprctile(cell_mass,cell_mass(cell_id)))]);
      cell_stage_cell_mass = eval(char(strcat('cell_mass_', cell_stage))); % within cell stage only
      text(0,.3,['Cell Size Percentile (within stage): ' num2str(invprctile(cell_stage_cell_mass,cell_mass(cell_id)))]);
      text(0,.2,['Average STRADa Percentile (Overall): ' num2str(invprctile(STRADa_total,STRADa_total(cell_id)))]);
      cell_stage_STRADa_total = eval(char(strcat('STRADa_total_', cell_stage))); % within cell stage only
      text(0,.1,['Average STRADa Percentile (within stage): ' num2str(invprctile(cell_stage_STRADa_total,STRADa_total(cell_id)))]);
      axis off

      %% Image of line scan from centroid to furthest point on cell boundry
      % get X any Y of each point on the perim
      [perim_y,perim_x]=find(bwperim(cyto_mask));
      % align the point so that it is relative to the nucleus's center, ie. it's origin (0,0) is the nucleus's center
      perim_x = perim_x - centroid_x;
      perim_y = perim_y - centroid_y;
      % get euclidean distance from manhatten distance
      [theta,rho] = cart2pol(perim_x,perim_y);
      [max_dist, max_dist_id] = max(rho);
      max_dist_x = perim_x(max_dist_id) + centroid_x;
      max_dist_y = perim_y(max_dist_id) + centroid_y;
      subplot(3,3,3)
      X = [centroid_x max_dist_x];
      Y = [centroid_y max_dist_y];
      line(X,Y,'Color','r')


      %% Plot line scan
      subplot(3,3,6)
      [CX,CY,C,xi,yi]=improfile(O.IM{cell_mass_channel},X,Y); 
      hold on
      plot(1:length(C),C,'-'); 
      grid on
      axis tight
      xlabel('lateral distance [px]'); 
      ylabel('intensity [a.u.]'); 
      title('STRADa linescan');

      % %% Plot line scan
      % subplot(3,3,9)
      % thickness=3;
      % [CX,CY,C_sum,C,xi,yi]=improfile_integrated(O.IM{cell_mass_channel},thickness,X,Y); 
      % if C_sum(1)<C_sum(end)
      %   C_sum=fliplr(C_sum);
      % end
      % hold on
      % plot(1:length(C),C,'.'); 
      % plot(1:length(C_sum),C_sum/(thickness+1),'-'); 
      % grid on
      % axis tight
      % xlabel('lateral distance [px]'); 
      % ylabel('intensity [a.u.]'); 
      % title('comparision of single linescan and integrated line scan');

      full_path = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 6 - normalized by total STRADa instead of SE/single cell images/ratio outliers/%s',aspects(i).foldername);
      mkdir(full_path)
      filename = sprintf('%s/ratio-%s_stage-%s_mass-%s_id-%s.png',full_path,num2str(STRADa_ratio(cell_id)),char(cell_stage),num2str(cell_mass(cell_id)),num2str(cell_id))
      export_fig(filename, '-m1')
  

      % Keep track of number saved images and quit if reached maximum
      if mod(c,2) == 0 %number is even
        high_count = high_count + 1;
      else
        low_count = low_count + 1;
      end
      
      if low_count >= max_count & high_count >= max_count
        l = 'done aspect'
        break
      end

    end
    


  end
end