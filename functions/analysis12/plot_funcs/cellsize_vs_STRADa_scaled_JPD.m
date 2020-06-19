function fun(ResultTable, Data, resolution)
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
  figure
  h=plot(X,Y,'.')

  ylabel('STRADa / Cell Size')
  ylabel(y_name,'Interpreter','none')
  set(h,'markersize',.01);

  hold on
  % Choose X and Y data to plot
  x_name = 'cell_size';
  y_name = 'STRADa_nuc ./ cell_size';

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

  % Plot
  h=plot(X,Y,'.')
  set(h,'markersize',.01);

  legend({'STRADa in Cytoplasm','STRADa in Nucleus'})
  legendmarkeradjust(30);

  % Save the plot to disk
  filename = sprintf('scaled_cyto_AND_nuc_STRADa_vs_cell_size.png');
  export_fig(['plots\' filename], '-m2')
end