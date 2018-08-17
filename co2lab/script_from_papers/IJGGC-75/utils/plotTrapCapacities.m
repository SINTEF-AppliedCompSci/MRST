function plotTrapCapacities(Gt, strap, p, t)
% Plots structural trapping capacity field (strap)
% Also shades region where CO2 is non-dense according to (p,t)


    % Map of structural traps and their capacities
    figure, title([num2str(sum(strap)/giga), ' Mt'])
    plotFaces(Gt, boundaryFaces(Gt));
    plotCellData(Gt, strap/giga, strap > 0, 'edgecolor','none');
    colorbar;
    axis equal tight off;
    set(gca,'FontSize',14);
    clim = get(gca,'CLim');

    
    % Include red shading to indicate where CO2 is non-dense as well as
    % non-supercritical
    hold on;
    suitable_cells = aquiferConditionSuitability(p, t, 'plot',false);
    non_suitable_cells = ~suitable_cells; 
    plotCellData(Gt, Gt.cells.z, non_suitable_cells, ...
        'facecolor','red', 'facealpha',0.2, 'edgecolor','none');
    set(gca,'CLim',clim);
    H = removeCells(Gt.parent, find(~non_suitable_cells));
    H = topSurfaceGrid(H);
    plotFaces(H, boundaryFaces(H), 'LineWidth',2, ...
                                   'FaceColor','r', ...
                                   'EdgeColor','r');


end