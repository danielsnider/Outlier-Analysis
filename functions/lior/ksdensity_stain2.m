function prob1d_by_rows(ResultTable, PlateMapRow)
%PROB1D_BY_ROWS Summary of this function goes here
%   Detailed explanation goes here
    
    % Find cells in each condition.
    Ctrl1 = ResultTable.Row==2 & ResultTable.Column==PlateMapRow.column; 
    Ctrl2 = ResultTable.Row==3 & ResultTable.Column==PlateMapRow.column;
    
    Low1 = ResultTable.Row==4 & ResultTable.Column==PlateMapRow.column;
    Low2 = ResultTable.Row==5 & ResultTable.Column==PlateMapRow.column;
    
    High1 = ResultTable.Row==6 & ResultTable.Column==PlateMapRow.column;
    High2 = ResultTable.Row==PlateMapRow.column & ResultTable.Column==PlateMapRow.column;
    
    % Convert the experiment row order from one big string to an array
    exp_row_order_array = strsplit(char(PlateMapRow.ExpRowOrder),',');

    % Plot each experiment in a loop
    figure;
    hold on;
    for exp_num=min(ResultTable.Row):max(ResultTable.Row)
        Cells = ResultTable.Row==exp_num & ResultTable.Column==PlateMapRow.column;
        exp_color = experiment_to_color(char(exp_row_order_array(exp_num)));
       [fx,sx]=ksdensity(ResultTable.NInt(Cells,PlateMapRow.stain2_channel_number), linspace(prctile(ResultTable.NInt(:, 2), .5), prctile(ResultTable.NInt(:,2), 95.5), 1000));
       plot(sx,fx,exp_color);
    end
    
    % Make the plot pretty
    xlabel(PlateMapRow.stain2_name);ylabel('Frequency');
    legend(exp_row_order_array(2:end));
    xlim([prctile(ResultTable.NInt(:, 2), .5) prctile(ResultTable.NInt(:,2), 95.5)]);
    
    % Save the plot to disk
    filename = sprintf('ksdensity_stain2_p%dc%d_%s.png',PlateMapRow.plate, PlateMapRow.column, char(PlateMapRow.stain2_name)); % example result p1c8_p21_pS6_JPD.png
    saveas(gcf,['plots\' filename]);
end
