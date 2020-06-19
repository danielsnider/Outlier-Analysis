function fun(ResultTable, Data, resolution)
  % Access data
  channels = Data.O.Channels;
  STRADa_channel = find(ismember(channels, 'STRADa'));
  cell_mass_channel = find(ismember(channels, 'SE'));
  DAPI_channel = find(ismember(channels, 'DAPI'));

  cell_size = ResultTable.CArea(:);
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;
  STRADa_nuc_scaled_for_size = STRADa_nuc./cell_size;
  STRADa_cyto_scaled_for_size = STRADa_cyto./cell_size;

  % Eliminate outliers
  lower_prctile = 1;
  higher_prctile = 99;
  cell_size = set_outliers_to(cell_size, lower_prctile, higher_prctile, NaN);
  STRADa_nuc = set_outliers_to(STRADa_nuc, lower_prctile, higher_prctile, NaN);
  STRADa_total = set_outliers_to(STRADa_total, lower_prctile, higher_prctile, NaN);
  STRADa_cyto = set_outliers_to(STRADa_cyto, lower_prctile, higher_prctile, NaN);
  STRADa_nuc_scaled_for_size = set_outliers_to(STRADa_nuc_scaled_for_size, lower_prctile, higher_prctile, NaN);
  STRADa_cyto_scaled_for_size = set_outliers_to(STRADa_cyto_scaled_for_size, lower_prctile, higher_prctile, NaN);


  % Choose X and Y data to plot
  x_name = 'cell_size';
  y_name = 'STRADa_cyto_scaled_for_size';
  z_name = 'STRADa_nuc_scaled_for_size';

  X=eval(x_name);
  Y=eval(y_name);
  Z=eval(z_name);

  % Eliminate some data
  X = X(1:100:end);
  Y = Y(1:100:end);
  Z = Z(1:100:end);

  gridx1 = linspace(prctile(X,1),prctile(X,99),resolution);
  gridx2 = linspace(prctile(Y,1),prctile(Y,99),resolution);
  gridx3 = linspace(prctile(Z,1),prctile(Z,99),resolution);
  [x1,x2,x3] = ndgrid(gridx1,gridx2,gridx3);
  x1 = x1(:,:)';
  x2 = x2(:,:)';
  x3 = x3(:,:)';
  xi = [x1(:) x2(:) x3(:)];
  three_dimensions = [ X, Y, Z ];
  [f,bw] = mvksdensity(three_dimensions,xi,'kernel','normpdf','Bandwidth',2);
  JPD = reshape(f,[resolution resolution resolution]);
  figure;
  imshow3Dfull(JPD,[]);
  colormap(gca,jet);
  % Make plot pretty
  
  set(gcf, 'XData', [min(X), max(X)], 'YData', [0, 100]);

  
  title('Nuc / Cyto')

  % Plot data
 % figure
 % h=scatter3(X,Y,Z,0.1,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

  xlabel(x_name)
  ylabel(y_name)
  zlabel(z_name)

  % Save the plot to disk
  filename = sprintf('cellsize_vs_STRADa_nuc_vs_STRADa_cyto_scatter3.png');
  export_fig(['plots\' filename], '-m2')
end