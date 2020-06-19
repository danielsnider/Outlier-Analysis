function cc_cond_by_rows(ResultTable, PlateMapRow)
%SCATTER_BY_ROWS Summary of this function goes here
%   Detailed explanation goes here


% % delete this after
% addpath '\\carbon.research.sickkids.ca\rkafri\Miriam\Matlab function library'
% addpath 'functions'
%
% PlateMap = readtable('Plate_Map.xlsx')
% % until here
% PlateMapRow = PlateMap(5,:)
% load(char(PlateMapRow.ResultsTableFilename))

% Return if there's not enough data for this plot
if strcmp(PlateMapRow.stain2_name,'') || strcmp(PlateMapRow.stain4_name,'')
    return
end

% Find cells in each condition.
Ctrl1 = ResultTable.Row==2 & ResultTable.Column==PlateMapRow.column;
Ctrl2 = ResultTable.Row==3 & ResultTable.Column==PlateMapRow.column;
Ctrl = (ResultTable.Row==3 | ResultTable.Row==2) & ResultTable.Column==PlateMapRow.column;

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
    
    %Cells = ResultTable.Row==exp_num & ResultTable.Column==PlateMapRow.column;
    
    
    %% CELL CYCLE CONDITIONAL PROBABILITY
    CellsCtrl1 = [ ResultTable.NInt(Ctrl1,2) ResultTable.CInt(Ctrl1,4)];
    CellsCtrl2 = [ ResultTable.NInt(Ctrl2,2) ResultTable.CInt(Ctrl2,4)];
    CellsLow1 = [ ResultTable.NInt(Low1,2) ResultTable.CInt(Low1,4)];
    CellsLow2 = [ ResultTable.NInt(Low2,2) ResultTable.CInt(Low2,4)];
    CellsHigh1 = [ ResultTable.NInt(High1,2) ResultTable.CInt(High1,4)];
    CellsHigh2 = [ ResultTable.NInt(High2,2) ResultTable.CInt(High2,4)];
    
    % Identify cell cycle stage based on DNA level
    [C1_nDNA, C1_idxG1, C1_idxS, C1_idxG2]=SortDNA(ResultTable.NInt(Ctrl1,1));
    [C2_nDNA, C2_idxG1, C2_idxS, C2_idxG2]=SortDNA(ResultTable.NInt(Ctrl2,1));
    [L1_nDNA, L1_idxG1, L1_idxS, L1_idxG2]=SortDNA(ResultTable.NInt(Low1,1));
    [L2_nDNA, L2_idxG1, L2_idxS, L2_idxG2]=SortDNA(ResultTable.NInt(Low2,1));
    [H1_nDNA, H1_idxG1, H1_idxS, H1_idxG2]=SortDNA(ResultTable.NInt(High1,1));
    [H2_nDNA, H2_idxG1, H2_idxS, H2_idxG2]=SortDNA(ResultTable.NInt(High2,1));
    
    
    % Make Meshgrid
    X = linspace(prctile(ResultTable.NInt(Ctrl,2),1),prctile(ResultTable.NInt(Ctrl,2),99),45);
    Y = linspace(prctile(ResultTable.CInt(Ctrl,4),1),prctile(ResultTable.CInt(Ctrl,4),99),45);
    [a, y] = meshgrid(X,Y);
    
    %Normalized DNA
    [c1_g1,sx]=ksdensity(ResultTable.NInt(Ctrl1(C1_idxG1)), linspace(prctile(ResultTable.NInt(Ctrl1(C1_idxG1)),1),prctile(ResultTable.NInt(Ctrl1(C1_idxG1)),99),2025));
    [c1_s,sx]=ksdensity(ResultTable.NInt(Ctrl1(C1_idxS)), linspace(prctile(ResultTable.NInt(Ctrl1(C1_idxS)),1),prctile(ResultTable.NInt(Ctrl1(C1_idxS)),99),2025));
    [c1_g2,sx]=ksdensity(ResultTable.NInt(Ctrl1(C1_idxG2)), linspace(prctile(ResultTable.NInt(Ctrl1(C1_idxG2)),1),prctile(ResultTable.NInt(Ctrl1(C1_idxG2)),99),2025));
    % PROBLEM ARISES HERE -------------------------------- the ResultsTable
    % is empty and messes up the ksdensity
    [c2_g1,sx]=ksdensity(ResultTable.NInt(Ctrl2(C2_idxG1)), linspace(prctile(ResultTable.NInt(Ctrl2(C2_idxG1)),1),prctile(ResultTable.NInt(Ctrl2(C2_idxG1)),99),2025));
    [c2_s,sx]=ksdensity(ResultTable.NInt(Ctrl2(C2_idxS)), linspace(prctile(ResultTable.NInt(Ctrl2(C2_idxS)),1),prctile(ResultTable.NInt(Ctrl2(C2_idxS)),99),2025));
    [c2_g2,sx]=ksdensity(ResultTable.NInt(Ctrl2(C2_idxG2)), linspace(prctile(ResultTable.NInt(Ctrl2(C2_idxG2)),1),prctile(ResultTable.NInt(Ctrl2(C2_idxG2)),99),2025));
    [l1_g1,sx]=ksdensity(ResultTable.NInt(Low1(L1_idxG1)), linspace(prctile(ResultTable.NInt(Low1(L1_idxG1)),1),prctile(ResultTable.NInt(Low1(L1_idxG1)),99),2025));
    [l1_s,sx]=ksdensity(ResultTable.NInt(Low1(L1_idxS)), linspace(prctile(ResultTable.NInt(Low1(L1_idxS)),1),prctile(ResultTable.NInt(Low1(L1_idxS)),99),2025));
    [l1_g2,sx]=ksdensity(ResultTable.NInt(Low1(L1_idxG2)), linspace(prctile(ResultTable.NInt(Low1(L1_idxG2)),1),prctile(ResultTable.NInt(Low1(L1_idxG2)),99),2025));
    [l2_g1,sx]=ksdensity(ResultTable.NInt(Low2(L2_idxG1)), linspace(prctile(ResultTable.NInt(Low2(L2_idxG1)),1),prctile(ResultTable.NInt(Low2(L2_idxG1)),99),2025));
    [l2_s,sx]=ksdensity(ResultTable.NInt(Low2(L2_idxS)), linspace(prctile(ResultTable.NInt(Low2(L2_idxS)),1),prctile(ResultTable.NInt(Low2(L2_idxS)),99),2025));
    [l2_g2,sx]=ksdensity(ResultTable.NInt(Low2(L2_idxG2)), linspace(prctile(ResultTable.NInt(Low2(L2_idxG2)),1),prctile(ResultTable.NInt(Low2(L2_idxG2)),99),2025));
    [h1_g1,sx]=ksdensity(ResultTable.NInt(High1(H1_idxG1)), linspace(prctile(ResultTable.NInt(High1(H1_idxG1)),1),prctile(ResultTable.NInt(High1(H1_idxG1)),99),2025));
    [h1_s,sx]=ksdensity(ResultTable.NInt(High1(H1_idxS)), linspace(prctile(ResultTable.NInt(High1(H1_idxS)),1),prctile(ResultTable.NInt(High1(H1_idxS)),99),2025));
    [h1_g2,sx]=ksdensity(ResultTable.NInt(High1(H1_idxG2)), linspace(prctile(ResultTable.NInt(High1(H1_idxG2)),1),prctile(ResultTable.NInt(High1(H1_idxG2)),99),2025));
    [h2_g1,sx]=ksdensity(ResultTable.NInt(High2(H2_idxG1)), linspace(prctile(ResultTable.NInt(High2(H2_idxG1)),1),prctile(ResultTable.NInt(High2(H2_idxG1)),99),2025));
    [h2_s,sx]=ksdensity(ResultTable.NInt(High2(H2_idxS)), linspace(prctile(ResultTable.NInt(High2(H2_idxS)),1),prctile(ResultTable.NInt(High2(H2_idxS)),99),2025));
    [h2_g2,sx]=ksdensity(ResultTable.NInt(High2(H2_idxG2)), linspace(prctile(ResultTable.NInt(High2(H2_idxG2)),1),prctile(ResultTable.NInt(High2(H2_idxG2)),99),2025));
    
    c1_g1_reshaped = reshape(c1_g1,[45 45]);
    c1_s_reshaped = reshape(c1_s,[45 45]);
    c1_g2_reshaped = reshape(c1_g2,[45 45]);
    c2_g1_reshaped = reshape(c2_g1,[45 45]);
    c2_s_reshaped = reshape(c2_s,[45 45]);
    c2_g2_reshaped = reshape(c2_g2,[45 45]);
    l1_g1_reshaped = reshape(l1_g1,[45 45]);
    l1_s_reshaped = reshape(l1_s,[45 45]);
    l1_g2_reshaped = reshape(l1_g2,[45 45]);
    l2_g1_reshaped = reshape(l2_g1,[45 45]);
    l2_s_reshaped = reshape(l2_s,[45 45]);
    l2_g2_reshaped = reshape(l2_g2,[45 45]);
    h1_g1_reshaped = reshape(h1_g1,[45 45]);
    h1_s_reshaped = reshape(h1_s,[45 45]);
    h1_g2_reshaped = reshape(h1_g2,[45 45]);
    h2_g1_reshaped = reshape(h2_g1,[45 45]);
    h2_s_reshaped = reshape(h2_s,[45 45]);
    h2_g2_reshaped = reshape(h2_g2,[45 45]);
    
    subplot(3,6,1) % control 1 conditionals given cell size stage
    [d, f]= ksdensity(CellsCtrl1, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(c1_g1_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('CONTROL 1: f%s,%s| G1)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    subplot(3,6,7) % control 1 conditionals given cell size stage
    [d, f]= ksdensity(CellsCtrl1, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(c1_s_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('CONTROL 1: f%s,%s| S)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    subplot(3,6,13) % control 1 conditionals given cell size stage
    [d, f]= ksdensity(CellsCtrl1, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(c1_g2_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('CONTROL 1: f%s,%s| G2)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    
    subplot(3,6,2) % controls 2 given cell size stage
    [d, f]= ksdensity(CellsCtrl2, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(c2_g1_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('CONTROL 2: f%s,%s| G1)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    subplot(3,6,8) % controls 2 given cell size stage
    [d, f]= ksdensity(CellsCtrl2, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(c2_s_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('CONTROL 2: f%s,%s| S)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    subplot(3,6,14) % controls 2 given cell size stage
    [d, f]= ksdensity(CellsCtrl2, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(c2_g2_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('CONTROL 2: f%s,%s| G2)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    
    subplot(3,6,3) % low 1 given cell size stage
    [d, f]= ksdensity(CellsLow1, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(l1_g1_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('LOW 1: f%s,%s| G1)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    subplot(3,6,9) % low 1 given cell size stage
    [d, f]= ksdensity(CellsLow1, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(l1_s_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('LOW 1: f%s,%s| S)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    subplot(3,6,15) % low 1 given cell size stage
    [d, f]= ksdensity(CellsLow1, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(l1_g2_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('LOW 1: f%s,%s| G2)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    
    subplot(3,6,4) % low 2 given cell size stage
    [d, f]= ksdensity(CellsLow2, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(l2_g1_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('LOW 2: f%s,%s| G1)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    subplot(3,6,10) % low 2 given cell size stage
    [d, f]= ksdensity(CellsLow2, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(l2_s_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('LOW 2: f%s,%s| S)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    subplot(3,6,16) % low 2 given cell size stage
    [d, f]= ksdensity(CellsLow2, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(l2_g2_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('LOW 2: f%s,%s| G2)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    
    subplot(3,6,5) % high 1 given cell size stage
    [d, f]= ksdensity(CellsHigh1, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(h1_g1_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('HIGH 1: f%s,%s| G1)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    subplot(3,6,11) % high 1 given cell size stage
    [d, f]= ksdensity(CellsHigh1, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(h1_s_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('HIGH 1: f%s,%s| S)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    subplot(3,6,17) % high 1 given cell size stage
    [d, f]= ksdensity(CellsHigh1, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(h1_g2_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('HIGH 1: f%s,%s| G2)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    
    subplot(3,6,6) % high 2 given cell size stage
    [d, f]= ksdensity(CellsHigh2, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(h2_g1_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('HIGH 2: f%s,%s| G1)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    subplot(3,6,12) % high 2 given cell size stage
    [d, f]= ksdensity(CellsHigh2, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(h2_s_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('HIGH 2: f%s,%s| S)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    subplot(3,6,18) % high 2 given cell size stage
    [d, f]= ksdensity(CellsHigh2, [a(:) y(:)]);
    d_reshaped = reshape(d,[45 45]);
    imagesc(d_reshaped ./repmat(sum(h2_g2_reshaped),45,1));colormap(gca,jet);set(gca,'ydir','normal')
    title('HIGH 2: f%s,%s| G2)',PlateMapRow.stain4_name,PlateMapRow.stain2_name );ylabel(PlateMapRow.stain4_name);xlabel(PlateMapRow.stain2_name);
    
end


% Save the plot to disk
filename = sprintf('cc_cond_p%dc%d_%s_%s.png',PlateMapRow.plate, PlateMapRow.column, char(PlateMapRow.stain2_name), char(PlateMapRow.stain4_name)); % example result p1c8_p21_pS6_JPD.png
saveas(gcf,['plots\' filename]);
end




