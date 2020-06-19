function prob1d_by_rows(ResultTable, PlateMapRow)
    
    % Convert the experiment row order from one big string to an array
    exp_row_order_array = strsplit(char(PlateMapRow.ExpRowOrder),',');
    if strcmp(PlateMapRow.stain2_name,'') || strcmp(PlateMapRow.stain4_name,'')
        return
    end 
    % Plot each experiment in a loop
    for exp_num=min(ResultTable.Row):max(ResultTable.Row)
        figure;
        Cells = ResultTable.Row==exp_num & ResultTable.Column==PlateMapRow.column;
        % Make a meshgrid for the heatmap
        X = linspace(prctile(ResultTable.NInt(Cells,2),1),prctile(ResultTable.NInt(Cells,2),99),200);
        Y = linspace(prctile(ResultTable.CInt(Cells,4),1),prctile(ResultTable.CInt(Cells,4),99),200);
        [a, y] = meshgrid(X,Y);
        exp_color = experiment_to_color(char(exp_row_order_array(exp_num)));
        Cells2 = [ ResultTable.NInt(Cells,2) ResultTable.CInt(Cells,4)];
        [d, f]= ksdensity(Cells2, [a(:) y(:)]);
        d_reshaped = reshape(d,[200 200]);
        imagesc(d_reshaped);colormap(gca,jet)
        set(gca,'ydir','normal')
    
        % Make the plot pretty
        xlabel(PlateMapRow.stain2_name);ylabel('pS6');
      %  legend(exp_row_order_array(2:end));
        %xlim([prctile(ResultTable.NInt(Cells, 2), .5) prctile(ResultTable.NInt(Cells,2), 95.5)]);
        % Save the plot to disk
        filename = sprintf('JPD_p%dc%d_%s_%s.png',PlateMapRow.plate, PlateMapRow.column, char(PlateMapRow.stain2_name), char(exp_row_order_array(exp_num))); % example result p1c8_p21_pS6_JPD.png
        saveas(gcf,['plots\' filename]);
    end
end
