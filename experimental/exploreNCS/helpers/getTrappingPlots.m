function [ capacityOutput, hfig, hax ] = getTrappingPlots(Gt, ta, rock2D, seaName)
% Structural trapping plots along with capacityOutput

% NB: put computations into helper functions to make code cleaner and
% easier to read.

    
    % Initialize figure:
    figure; set(gcf,'Position',[1 1 3000 500])
    hfig = gcf;
    
    
    % ---------------------------------------------------------------------
    % PLOT TRAPS COLORED BY CO2 MASS STORAGE CAPACITY
    
    % bulk (rock + pore) volume of trap
    bulkVols = volumesOfTraps(Gt, ta);
    
    % pore volume of trap
    numTraps = max(unique(ta.traps));
    assert(numel(bulkVols)==numTraps,'Mis-match of total trap number.');
    trap_ixs = 1:numTraps;
    poreVols = computeTrapVolume(Gt, ta, rock2D.poro, trap_ixs);
    assert(numel(bulkVols) == numel(poreVols),'Mis-match of bulk and pore volumes.')
    
    
    % Some output:
    fprintf('  Num. global traps: %d\n', numel(bulkVols));
    fprintf('  Total trap (bulk) volume: %e m3\n', sum(bulkVols));
    fprintf('  Avg. trap (bulk) size: %e m3\n', mean(bulkVols));
    fprintf('  Total trap pore volume: %e m3\n', sum(poreVols));
    fprintf('  Avg. trap pore size: %e m3\n', mean(poreVols));
    
    
    % TRAP VOLUME
    subplot(1,6,1); hsub1 = gca; hfsub1 = gcf;
    hold on
    bf = boundaryFaces(Gt);
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    trapcells = ta.traps~=0;
    cellsTrapVol = zeros(Gt.cells.num,1);
    cellsTrapVol(trapcells) = bulkVols(ta.traps(trapcells));
    plotCellData(Gt, cellsTrapVol/1e3/1e3/1e3, cellsTrapVol~=0, 'EdgeColor','none')
    title({'Trap Volume';'(km^3)'}); axis off; colorbar
    
   
    % TRAP PORE VOLUME
    subplot(1,6,2); hsub2 = gca; hfsub2 = gcf;
    hold on
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    cellsTrapPoreVol = zeros(Gt.cells.num,1);
    cellsTrapPoreVol(trapcells) = poreVols(ta.traps(trapcells));
    plotCellData(Gt, cellsTrapPoreVol/1e3/1e3/1e3, cellsTrapPoreVol~=0, 'EdgeColor','none')
    title({'Trap Pore Volume';'(km^3)'}); axis off; colorbar
    
    
    % ---------------------------------------------------------------------
    % GET TRAPPING BREAKDOWN: structural, residual, dissolution
    % first, compute theoretical capacity (upper bound):
    info = getSeaInfo(seaName);
    rhoCref = info.water_density;
    [ capacityOutput ] = getTrappingCapacities(Gt, rock2D, ta, ...
        rhoCref, info.water_density, info.seafloor_temp, info.seafloor_depth, ...
        info.temp_gradient, info.press_deviation, info.res_sat_co2, info.res_sat_wat, ...
        info.dis_max);

    % Distributed CO2 mass under structural traps: 
    cellsTrapCO2Mass = zeros(Gt.cells.num,1);
    cellsTrapCO2Mass(trapcells) = capacityOutput.strap_mass_co2(trapcells);

    % Cumulative CO2 mass under structural traps:
    trapcaps = accumarray(ta.traps(trapcells), capacityOutput.strap_mass_co2(trapcells));
    trapcap_tot = zeros(Gt.cells.num,1); %ones(size(ta.traps)) * NaN;
    trapcap_tot(trapcells) = trapcaps(ta.traps(trapcells));
    
    
    % DISTRIBUTED CO2 MASS
    subplot(1,6,3); hsub3 = gca; hfsub3 = gcf;
    hold on
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt, cellsTrapCO2Mass/1e9, cellsTrapCO2Mass~=0, 'EdgeColor','none')
    title({'Distributed CO2 Mass';'under Trap (Mt)'}); axis off; colorbar;

    
    % ACCUMULATED CO2 MASS
    subplot(1,6,4); hsub4 = gca; hfsub4 = gcf;
    hold on
    %plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'none')
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt, trapcap_tot/1e9, trapcap_tot~=0, 'EdgeColor','none')
    title({'Accumulated CO2 Mass';'under Trap (Mt)'}); axis off; colorbar;


    % ---------------------------------------------------------------------
    % PLOT REACHABLE CAPACITY
    trees = maximizeTrapping(Gt, 'res', ta, 'calculateAll', true, 'removeOverlap', false);
    tvols = [trees.value]; %#ok
    int_tr = find(ta.trap_regions); %#ok ixs of cells spilling into interior trap
    [dummy, reindex] = sort([trees.root], 'ascend'); %#ok
    structural_mass_reached = zeros(Gt.cells.num, 1);
    for i = 1:numel(ta.trap_z) % loop over each trap
        % ix of cells spilling directly into this trap
        cix = find(ta.trap_regions == i);
        % cell indices of all cells of this trap, and its upstream traps
        aix = find(ismember(ta.traps, [trees(reindex(i)).traps]));
        % compute total structural trap capacity (in mass terms) of these
        % cells, and store result
        structural_mass_reached(cix) = sum(capacityOutput.strap_mass_co2(aix)); %#ok
    end

    %
    subplot(1,6,5); hsub5 = gca; hfsub5 = gcf;
    hold on
    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt, structural_mass_reached/1e3/1e6, structural_mass_reached~=0, 'EdgeColor','none');
    title({'Reachable structural';'capacity (Mt)'}); axis off; colorbar;


    % ---------------------------------------------------------------------
    % PLOT SPILL PATHS AND TOPOLOGY
    subplot(1,6,6); hsub6 = gca; hfsub6 = gcf;
    hold on
    mapPlot(gcf, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines);
    title({'Topology, traps';'and spill-paths'}); axis off
    

    
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    % Adjust axis, fonts
    hfig = gcf;
    set(findobj(hfig.Children,'Type','axes'),'Fontsize',16);
    set(findobj(hfig.Children,'Type','colorbar'),'Fontsize',16);
    axis(findobj(hfig.Children,'Type','axes'),'equal','tight');
    
    % After tighting of axis, re-size mapPlot since the absence of a
    % colorbar causes subplot window to be larger than other subplots
    hfig.Children(1).Position(2) = hfig.Children(3).Position(2);
    hfig.Children(1).Position(3) = hfig.Children(3).Position(3);
    hfig.Children(1).Position(4) = hfig.Children(3).Position(4);

    
    hfig = gcf;
    hax  = gca;
    
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    % Capacity Output
    capacityOutput.structural_mass_reached = structural_mass_reached/1e3/1e6; % Mt


end