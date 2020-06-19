function Show_Selected_Cell(ResultTable,Data)
  %% Plot data, choose cell, and display the image and boundry of that cell.

  % Access data
  channels = Data.O.Channels;
  STRADa_channel = find(ismember(channels, 'STRADa'));
  cell_mass_channel = find(ismember(channels, 'SE'));
  DAPI_channel = find(ismember(channels, 'DAPI'));
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;
  cell_size = ResultTable.CArea(:);

  % Choose X and Y data to plot
  x_name = 'cell_size';
  y_name = 'STRADa_cyto ./ cell_size';

  X=eval(x_name);
  Y=eval(y_name);


  % Eliminate outliers
  lower_prctile=1;
  higher_prctile=99;
  X_lower_prctile = prctile(X,lower_prctile);
  X_higher_prctile = prctile(X,higher_prctile);
  Y_lower_prctile = prctile(Y,lower_prctile);
  Y_higher_prctile = prctile(Y,higher_prctile);
  X(X < X_lower_prctile | X > X_higher_prctile) = NaN;
  Y(Y < Y_lower_prctile | Y > Y_higher_prctile) = NaN;

  % Plot data
  figure(136)
  plot(X,Y,'.')

  xlabel(x_name,'Interpreter','none')
  ylabel(y_name,'Interpreter','none')

  % Choose cell of interest by clicking data point in plot
  [x,y,button] = ginput(1);

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