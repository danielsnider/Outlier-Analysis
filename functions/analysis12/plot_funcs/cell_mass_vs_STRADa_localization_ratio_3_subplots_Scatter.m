function fun(ResultTable, Data, resolution)
  channels = Data.O.Channels;
  cell_size_channel = find(ismember(channels, 'SE'));
  STRADa_channel = find(ismember(channels, 'STRADa'));
  cell_size = ResultTable.CInt(:,cell_size_channel);
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;
  STRADa_localization_ratio_nuc_cyto = STRADa_nuc./STRADa_cyto;
  STRADa_localization_ratio_nuc_total = STRADa_nuc./STRADa_total;
  STRADa_localization_ratio_cyto_total = STRADa_cyto./STRADa_total;

  all_STRA = [STRADa_localization_ratio_nuc_cyto ...
              STRADa_localization_ratio_nuc_total ...
              STRADa_localization_ratio_cyto_total]; % used for consistent plot axis limits

  % Scatter plot
  figure('Position',[251,    91,   1583,   708]);
  
  %% STRADa_nuc./STRADa_cyto
  subplot(1,3,1)
  idx = isfinite(cell_size) & isfinite(STRADa_localization_ratio_nuc_cyto);
  fit_line = fit(cell_size(idx),STRADa_localization_ratio_nuc_cyto(idx),'poly1');
  h=plot(fit_line,cell_size, STRADa_localization_ratio_nuc_cyto, '.b');
  % Make plot pretty
  title('Nuc / Cyto')
  legend(gca,'off');
  set(h,'markersize',.01);
  xlabel('Cell Size (SE)');ylabel('STRADa_n_u_c/STRADa_c_y_t_o');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(all_STRA(:), 1) prctile(all_STRA(:), 99)]);

  
  %% STRADa_nuc./STRADa_total
  subplot(1,3,2)
  idx = isfinite(cell_size) & isfinite(STRADa_localization_ratio_nuc_total);
  fit_line = fit(cell_size(idx),STRADa_localization_ratio_nuc_total(idx),'poly1');
  h=plot(fit_line,cell_size, STRADa_localization_ratio_nuc_total, '.b');
  % Make plot pretty
  title('Nuc / Total')
  legend(gca,'off');
  set(h,'markersize',.01);
  xlabel('Cell Size (SE)');ylabel('STRADa_n_u_c/STRADa_t_o_t_a_l');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(all_STRA(:), 1) prctile(all_STRA(:), 99)]);

  
  %% STRADa_cyto./STRADa_total
  subplot(1,3,3)
  idx = isfinite(cell_size) & isfinite(STRADa_localization_ratio_cyto_total);
  fit_line = fit(cell_size(idx),STRADa_localization_ratio_cyto_total(idx),'poly1');
  h=plot(fit_line,cell_size, STRADa_localization_ratio_cyto_total, '.b');
  % Make plot pretty
  title('Cyto / Total')
  legend(gca,'off');
  set(h,'markersize',.01);
  xlabel('Cell Size (SE)');ylabel('STRADa_c_y_t_o/STRADa_t_o_t_a_l');
  xlim([prctile(cell_size, 1) prctile(cell_size, 99)]);
  ylim([prctile(all_STRA(:), 1) prctile(all_STRA(:), 99)]);

  % Save the plot to disk
  filename = sprintf('3scatter.png');
  saveas(gcf,['plots\' filename]);
end

