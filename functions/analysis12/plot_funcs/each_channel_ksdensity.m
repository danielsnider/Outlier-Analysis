function fun(ResultTable, Data)
  channels = Data.O.Channels;
  dataset_name = Data.O.SegmentationParameters.DataSetName;

  % NUC
  figure('Position',[251,    91,   683,   708]);
  
  % Plot each channel in a loop
  for i=1:length(channels)
    subplot(2,2,i)
    [fx,sx]=ksdensity(ResultTable.NInt(:,i), linspace(prctile(ResultTable.NInt(:, i), .5), prctile(ResultTable.NInt(:, i), 95.5), 1000));
    plot(sx,fx);
    title(channels{i})
  end
  
  % Titles
  [ax1,h1]=suplabel('Single Cell Channel Intensity in Nucleus');
  set(h1,'FontSize',14)
  [ax2,h2]=suplabel('Frequency','y');
  set(h2,'FontSize',14)
  [ax3,h3]=suplabel(['Dataset: ' dataset_name],'t');
  set(h3,'FontSize',15)

  % Save the plot to disk
  filename = sprintf('ksdensity_each_channel_nuc.png');
  saveas(gcf,['plots\' filename]);

  %% CYTO
  figure('Position',[251,    91,   683,   708]);
  
  % Plot each channel in a loop
  for i=1:length(channels)
    subplot(2,2,i)
    total_intensity = ResultTable.CInt(:,i);
    nuc_intensity = ResultTable.NInt(:,i);
    cyto_intensity = total_intensity - nuc_intensity;
    [fx,sx]=ksdensity(cyto_intensity, linspace(prctile(cyto_intensity, .5), prctile(cyto_intensity, 95.5), 1000));
    plot(sx,fx);
    title(channels{i})
  end
  
  % Titles
  [ax1,h1]=suplabel('Single Cell Channel Intensity in Cytoplasm');
  set(h1,'FontSize',14)
  [ax2,h2]=suplabel('Frequency','y');
  set(h2,'FontSize',14)
  [ax3,h3]=suplabel(['Dataset: ' dataset_name],'t');
  set(h3,'FontSize',15)

  % Save the plot to disk
  filename = sprintf('ksdensity_each_channel_cyto.png');
  saveas(gcf,['plots\' filename]);

  %% log of STRADa in Nuc 
  figure('Position',[251,    91,   683,   708]);

  channels = Data.O.Channels;
  STRADa_channel = find(ismember(channels, 'STRADa'));
  STRADa_nuc_log = log(ResultTable.NInt(:,STRADa_channel));
  [fx,sx]=ksdensity(STRADa_nuc_log, linspace(prctile(STRADa_nuc_log, .5), prctile(STRADa_nuc_log, 95.5), 1000));
  plot(sx,fx);
  
  % Titles
  xlabel('Single Cell STRADa log(Intensity) in Nucleus');
  ylabel('Frequency');

  % Save the plot to disk
  filename = sprintf('ksdensity_STRADa_nuc_log.png');
  saveas(gcf,['plots\' filename]);
end
