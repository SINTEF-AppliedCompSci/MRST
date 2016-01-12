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
plotsOn = true;
figDirName = 'HammerfestProspectFigs';
mkdir(figDirName)

%% Get Sto dataset
fmName = 'Stofm';
[Gt, rock2D] = getFormationTopGrid(fmName, N);
bfinx        = getFormationsClosedBdryFaces(fmName, Gt);
rhoCref      = 760 * kilogram/meter^3;
seainfo      = getSeaInfo('BarentsSea', rhoCref);
trapInfo     = getTrappingInfo(Gt, rock2D, seainfo, 'plotsOn',false, ...
                        'fmName',fmName, 'closed_boundary_edges',bfinx);

% re-name variables
Gt_st       = Gt;
G_st        = Gt_st.parent;
ta_st       = trapInfo.ta;
rock2D_st   = rock2D;


%% Get cell index of Sto grid corresponding to NPD's prospects/closures
[ Grids, pts ] = getProspectTrapIDsAndGrids( G_st, Gt_st, ta_st );



%% Get volume and capacity of prospects, for comparison to NPD (Atlas chp 6)


names       = fieldnames(Grids);
capacities  = cell(1,numel(names));
for i = 1:numel(names)
    
    G           = getfield(Grids, names{i});
    trapcells   = G.trapcells;
    trapName    = { regexprep(G.name,'[^\w'']','') };  % cell for table
    
    [~, capacities{i}] = getTrappingInfo(Gt, rock2D, seainfo, 'plotsOn',false, ...
                                   'cells',trapcells, 'trapName',trapName, ...
                                   'closed_boundary_edges',bfinx);
                               
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
figure; set(gcf,'Position',[653 453 820 819])
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
% NB: annotation location is in normalized units for figure, while text is
% in the coordinate data units
annotation('textarrow',[0.509756097560976 0.569512195121951],...
    [0.743589743589744 0.663003663003663],'String','D','FontSize',16);
annotation('textarrow',[0.567073170731707 0.619512195121951],...
    [0.774114774114774 0.697191697191697],'String','C','FontSize',16);
annotation('textarrow',[0.707317073170732 0.601219512195122],...
    [0.477411477411477 0.527472527472527],'String','E','FontSize',16);
annotation('textarrow',[0.692682926829268 0.615853658536585],...
    [0.413919413919414 0.454212454212454],'String','F','FontSize',16);
annotation('textarrow',[0.48780487804878 0.44129737776251],...
    [0.258852258852259 0.285288270377731],'String','G','FontSize',16);
annotation('textarrow',[0.360975609756098 0.313832589030116],...
    [0.129426129426129 0.153167684976829],'String','H','FontSize',16);
annotation('textarrow',[0.318292682926829 0.41090976754838],...
    [0.63003663003663 0.552910862055991],...
    'String',['Greater',sprintf('\n'),'Sn\ohvit'],...
    'FontSize',16);
annotation('textarrow',[0.230487804878049 0.36570755102052],...
    [0.401709401709402 0.398561897567862],...
    'String',['Greater',sprintf('\n'),'Askeladd'],...
    'FontSize',16);
annotation('textarrow',[0.254878048780488 0.447365613136554],...
    [0.532356532356532 0.461960973640893],...
    'String',['Greater',sprintf('\n'),'Albatross'],...
    'FontSize',16);
hfig = gcf;
set(findall(hfig,'Type','TextArrow'), 'FontSize',18, 'LineWidth',1, ...
    'TextColor',[0 0 1], 'Color',[0 0 1])

if plotsOn
    pause
    export_fig(gcf,[figDirName '/' 'ProspectsInHammerfest_ref',num2str(N)], '-png','-transparent')
    %close
end

% text(pts.Xpt_gsf,  pts.Ypt_gsf,  'G. Sn\ohvit',   'HorizontalAlignment','center');
% text(pts.Xpt_alba, pts.Ypt_alba, 'G. Albatross', 'HorizontalAlignment','right');
% text(pts.Xpt_aske, pts.Ypt_aske, 'G. Askeladd',  'HorizontalAlignment','center');
% text(pts.Xpt_prosC, pts.Ypt_prosC,   'C', 'HorizontalAlignment','center');
% text(pts.Xpt_prosD, pts.Ypt_prosD,   'D', 'HorizontalAlignment','right');
% text(pts.Xpt_prosE, pts.Ypt_prosE,   'E', 'HorizontalAlignment','center');
% text(pts.Xpt_prosF, pts.Ypt_prosF,   'F', 'HorizontalAlignment','center');
% text(pts.Xpt_prosG, pts.Ypt_prosG,   'G', 'HorizontalAlignment','center');
% text(pts.Xpt_prosH, pts.Ypt_prosH,   'H', 'HorizontalAlignment','right');
% hax = gca;
% for i=1:numel(hax.Children)
%     if strcmpi(class(hax.Children(i)),'matlab.graphics.primitive.Text')
%         set(hax.Children(i), 'FontSize',18, 'FontWeight','bold', ...
%             'Color','blue', 'FontAngle','italic');
%     end
% end


% plot to show trapping capacity (Mt) of each prospect/closure
figure; set(gcf,'Position',[653 453 820 819])
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
if plotsOn
    export_fig(gcf,[figDirName '/' 'ProspectCapacity_ref',num2str(N)], '-png','-transparent')
    %close
end


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
for i = 1:numel(capacities)
    Seff = getStorageEfficiency(names{i}); % percent
    adjustedCap_Mt(i) = capacities{1,i}.storageCap_Mt * Seff/100;
end



% TODO:
% - convert this script into a function to include all helper functions
% within?
% - get trapping capacities in formation that contains faults



