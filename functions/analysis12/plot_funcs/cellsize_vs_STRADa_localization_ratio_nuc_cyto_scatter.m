function fun(ResultTable, Data, special)
  channels = Data.O.Channels;
  STRADa_channel = find(ismember(channels, 'STRADa'));
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
  lower_prctile = .1;
  higher_prctile = 99.9;
  cell_size = set_outliers_to(cell_size, lower_prctile, higher_prctile, NaN);
  STRADa_nuc = set_outliers_to(STRADa_nuc, lower_prctile, higher_prctile, NaN);
  STRADa_total = set_outliers_to(STRADa_total, lower_prctile, higher_prctile, NaN);
  STRADa_cyto = set_outliers_to(STRADa_cyto, lower_prctile, higher_prctile, NaN);
  STRADa_localization_ratio_nuc_cyto = set_outliers_to(STRADa_localization_ratio_nuc_cyto, lower_prctile, higher_prctile, NaN);
  STRADa_localization_ratio_nuc_total = set_outliers_to(STRADa_localization_ratio_nuc_total, lower_prctile, higher_prctile, NaN);
  STRADa_localization_ratio_cyto_total = set_outliers_to(STRADa_localization_ratio_cyto_total, lower_prctile, higher_prctile, NaN);

  % Choose X and Y data to plot
  x_name = 'cell_size';
  y_name = 'STRADa_localization_ratio_nuc_cyto';

  X=eval(x_name);
  Y=eval(y_name);


  figure('Position',[251,    91,   883,   708]);
  
  if strcmp(special, 'log')
    Y = log(Y);
  end
  h=plot(X,Y,'.')
  set(h,'markersize',.01);

  % Make plot pretty
  title('Nuc / Cyto')
  xlabel(x_name,'Interpreter','none');
  ylabel(y_name,'Interpreter','none');

  % Save the plot to disk
  filename = sprintf('plots/cellsize_vs_STRADa_localization_ratio_nuc_cyto_scatter_%s.png',special);
  export_fig(filename, '-m2')
end

