%% explore Prospects and Structural Closures in Hammerfest Basin Aquifer
% The NPD has tabulated storage capacity estimates for several regions that
% are found within the Hammerfest Basin aquifer, most notably in the Sto
% formation. These regions are called prospects or structural closures.
% Here, we extract the prospects and structural closures from the Sto
% formation, by assuming they correspond to the structural traps found by
% the MRST algorithm of trapAnalysis(). Two prospects are found within one
% of the structural closure, so to extract those prospects, an additional
% trapping analysis is performed on the structural closure.

% Note, prospect G and H are not quite captured by an entire structural
% trap, and may require other properties to extract them from the Sto
% formation grid (such as elevation or permeability properties).


mrstModule add libgeometry opm_gridprocessing

N = 1;

%% Get Sto dataset
Sto_trapInfo = getTrappingInfo('Stofm', N, 'plotsOn',false);

% re-name variables
Gt_st       = Sto_trapInfo.Gt;
G_st        = Gt_st.parent;
ta_st       = Sto_trapInfo.ta;
rock2D_st   = Sto_trapInfo.rock2D;


%% Get cell index of Sto grid corresponding to NPD's prospects/closures
[ Grids, pts ] = getProspectTrapIDsAndGrids( G_st, Gt_st, ta_st );



%% Get volume and capacity of prospects, for comparison to NPD (Atlas chp 6)
rhoCref = 760 * kilogram/meter^3;
info    = getSeaInfo('BarentsSea', rhoCref);

names       = fieldnames(Grids);
capacities  = cell(1,numel(names));
for i = 1:numel(names)
    
    G           = getfield(Grids, names{i});
    trapcells   = G.trapcells;
    trapName    = { regexprep(G.name,'[^\w'']','') };  % cell for table
    
    [~, capacities{i}] = getTrappingInfo('Stofm',1, 'plotsOn',false, ...
                                   'cells',trapcells, 'trapName',trapName);
                               
end
    
% Merge tables to get combined table of data:
Tcurrent        = capacities{1}.T;
% for latex table format
ColOfAmpSym     = table( repmat('&',numel(Tcurrent.Properties.RowNames),1), ...
    'RowNames',Tcurrent.Properties.RowNames, 'VariableNames',{'AmpSym'}); 
ColOfSlashSym   = table( repmat('\\',numel(Tcurrent.Properties.RowNames),1), ...
    'RowNames',Tcurrent.Properties.RowNames, 'VariableNames',{'SlashSym'}); 
Tcurrent        = join(ColOfAmpSym, Tcurrent, 'Keys','RowNames');
for i = 2:numel(capacities)
    T2add = capacities{i}.T;
    T2add = join( ColOfAmpSym, T2add, 'Keys','RowNames');
    Tcurrent = join( Tcurrent, T2add, 'Keys','RowNames');
end
Tcurrent = join( Tcurrent, ColOfSlashSym, 'Keys','RowNames');

% Display summary table in command window
Tcurrent %#ok


%% Colorize Prospects/Structural closures by storage capacity (Mt)

% plot to show all structrual traps identified by trapAnalysis, and then
% the boundaries of the prospect regions ontop to show comparison
figure
mapPlot(gcf, Gt_st, 'traps',ta_st.traps)
names = fieldnames(Grids);
for i = 1:numel(names)
    G = getfield(Grids, names{i});
    Gt = topSurfaceGrid(G);
    % plot boundaries of regions
    plotFaces(Gt, boundaryFaces(Gt), [1 0 0], 'Linewidth',3)
end
axis equal tight off
% add prospect labels on top of plot
text(pts.Xpt_gsf,  pts.Ypt_gsf,  'G. Sn\ohvit',   'HorizontalAlignment','center');
text(pts.Xpt_alba, pts.Ypt_alba, 'G. Albatross', 'HorizontalAlignment','right');
text(pts.Xpt_aske, pts.Ypt_aske, 'G. Askeladd',  'HorizontalAlignment','center');
text(pts.Xpt_prosC, pts.Ypt_prosC,   'C', 'HorizontalAlignment','center');
text(pts.Xpt_prosD, pts.Ypt_prosD,   'D', 'HorizontalAlignment','right');
text(pts.Xpt_prosE, pts.Ypt_prosE,   'E', 'HorizontalAlignment','center');
text(pts.Xpt_prosF, pts.Ypt_prosF,   'F', 'HorizontalAlignment','center');
text(pts.Xpt_prosG, pts.Ypt_prosG,   'G', 'HorizontalAlignment','center');
text(pts.Xpt_prosH, pts.Ypt_prosH,   'H', 'HorizontalAlignment','right');
hax = gca;
for i=1:numel(hax.Children)
    if strcmpi(class(hax.Children(i)),'matlab.graphics.primitive.Text')
        set(hax.Children(i), 'FontSize',18, 'FontWeight','bold', ...
            'Color','blue', 'FontAngle','italic');
    end
end


% plot to show trapping capacity (Mt) of each prospect/closure
figure;
plotCellData(Gt_st, ta_st.traps, 'EdgeColor','none','FaceColor',[1 0.9 0.9])
names = fieldnames(Grids);
for i = 1:numel(names)
    G = getfield(Grids, names{i});
    Gt = topSurfaceGrid(G);
    
    % colorize by storage capacity
    % TODO: ensure Grid and Capacity is for same trapName
    plotCellData(Gt, repmat(capacities{i}.storageCap_Mt, [Gt.cells.num,1]), 'EdgeColor','none')
    fprintf('\n Grid %s has storage capacity of %6.1f Mt. \n', names{i}, capacities{i}.storageCap_Mt);
    
    % plot boundaries of regions
    plotFaces(Gt, boundaryFaces(Gt), 'Linewidth',3)
end
light('Position',[-1 -1 1],'Style','infinite');lighting phong
axis equal tight off
hcb = colorbar; set(hcb,'FontSize',14);
set(hcb.Label, 'String','Mt', 'FontSize',16, 'Rotation',0);
set(gca,'DataAspect',[1 1 0.02]);
%view([-100 55]);
view(2)


% % using mapPlot to show contours --> TODO: prevent colors from being
% % changed by mapPlot
% figure;
% plotCellData(Gt_st, ta_st.traps, 'EdgeColor','none','FaceColor',[1 0.9 0.9])
% mapPlot(gcf, Gt_st);
% names = fieldnames(Grids);
% for i = 1:numel(names)
%     G = getfield(Grids, names{i});
%     Gt = topSurfaceGrid(G);
%     
%     % colorize by storage capacity
%     % TODO: ensure Grid and Capacity is for same trapName
%     plotCellData(Gt, repmat(capacities{i}.storageCap_Mt, [Gt.cells.num,1]), 'EdgeColor','none')
%     
%     % plot boundaries of regions
%     plotFaces(Gt, boundaryFaces(Gt), 'Linewidth',3)
% end
% mapPlot(gcf, Gt_st);
% %light('Position',[-1 -1 1],'Style','infinite');lighting phong





%% Use NPD's storage efficiency values to get reduced storage capacity estimates



% other things to try:
% - increase resolution of initial grid, to get traps again
% - compare between cell-based and node-based method of trapAnalysis
% - use lower CO2 density (700 kg/m3) as used by NPD calculations
% - use formula with constant CO2 density to compute storage capacity


