function prob1d_by_rows(ResultTable, PlateMapRow)
    
    % Convert the experiment row order from one big string to an array
    exp_row_order_array = strsplit(char(PlateMapRow.ExpRowOrder),',');
    if strcmp(PlateMapRow.stain2_name,'') || strcmp(PlateMapRow.stain4_name,'')
        return
    end 
    
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
%         imagesc(d_reshaped);colormap(gca,jet)
%         set(gca,'ydir','normal')
        
        xlabel('p-H2AX freq ');ylabel('PS6 freq');
        % conditional prob
        subplot (1,2,1);
        probMTORgivenDDR = d_reshaped ./repmat(sum(d_reshaped), 200, 1);
        imagesc(probMTORgivenDDR);colormap(gca,jet);
        title('CONTROL: f(mtor|ddr)') 
        set(gca,'ydir','normal')
        xlabel(PlateMapRow.stain2_name);ylabel('pS6');
        
        subplot (1,2,2);
        probDDRgivenMTOR = d_reshaped ./repmat(sum(d_reshaped'), 200, 1);
        imagesc(probDDRgivenMTOR);colormap(gca,jet);
        title('CONTROL: f(ddr|mtor)')
        set(gca,'ydir','normal')
        
        % Make the plot pretty
        xlabel(PlateMapRow.stain2_name);ylabel('pS6');
        
        % Save the plot to disk
        filename = sprintf('cond_p%dc%d_%s_%s.png',PlateMapRow.plate, PlateMapRow.column, char(PlateMapRow.stain2_name), char(exp_row_order_array(exp_num))); % example result p1c8_p21_pS6_JPD.png
        saveas(gcf,['plots\' filename]);
    end
end
