function contour2d_by_rows(ResultTable, PlateMapRow)
%SCATTER_BY_ROWS Summary of this function goes here
%   Detailed explanation goes here

    % Return if there's not enough data for this plot
    if strcmp(PlateMapRow.stain2_name,'') || strcmp(PlateMapRow.stain4_name,'')
        return
    end 

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
    for exp_num=2:7
        
        Cells = ResultTable.Row==exp_num & ResultTable.Column==PlateMapRow.column;
        X = linspace(prctile(ResultTable.NInt(Cells,2),1),prctile(ResultTable.NInt(Cells,2),99),200);
        Y = linspace(prctile(ResultTable.CInt(Cells,4),1),prctile(ResultTable.CInt(Cells,4),99),200);
        [a, y] = meshgrid(X,Y);
        
        exp_color = experiment_to_color(char(exp_row_order_array(exp_num)));
        Cells2 = [ ResultTable.NInt(Cells,2) ResultTable.CInt(Cells,4)];
        [d, f]= ksdensity(Cells2, [a(:) y(:)]);
        % SCATTER PLOT
        h=plot(ResultTable.NInt(Cells,PlateMapRow.stain2_channel_number),...
            ResultTable.CInt(Cells,PlateMapRow.stain4_channel_number),...
            ['.' exp_color]);
        set(h,'markersize',.01);
        
        % CONTOUR LINES 
        [d, f, bwc] = ksdensity(Cells2, [a(:) y(:)]);
        Z = reshape(d,[200 200]);
        contour(X,Y,Z,20,'LineColor',exp_color)
    end
    
    % Make the plot pretty
    xlabel(PlateMapRow.stain2_name);ylabel(PlateMapRow.stain4_name);
    legend(exp_row_order_array(2:7))
    xlim([prctile(ResultTable.NInt(:, 2), .5) prctile(ResultTable.NInt(:,2), 99.5)]);
    ylim([prctile(ResultTable.CInt(:,4), .5) prctile(ResultTable.CInt(:,4), 99.5)]);
    % Save the plot to disk
    filename = sprintf('contour2d_p%dc%d_%s_%s.png',PlateMapRow.plate, PlateMapRow.column, char(PlateMapRow.stain2_name), char(PlateMapRow.stain4_name)); % example result p1c8_p21_pS6_JPD.png
    saveas(gcf,['plots\' filename]);
end