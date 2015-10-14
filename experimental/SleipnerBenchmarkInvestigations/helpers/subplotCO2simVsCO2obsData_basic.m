function [ hfig, hax ] = subplotCO2simVsCO2obsData_basic(Years2plot, inj_year, plume, sim_report, Gt, states, fluid, model, wellXcoord, wellYcoord, wellCoord_x, wellCoord_y, trapstruct, ZoomIntoPlume, CO2plumeOutline_SatTol)

% Reservoir Time2plot is an array of all times in seconds you wish to
% visualize the CO2 plume (saturation). If any times coorespond to the
% observation plume data, that plume outline is plotted as well.

ReservoirTime2plot  = (Years2plot - inj_year(1)+1 ).*(365*24*60*60); % seconds

maxMassCO2 = zeros(1,numel(ReservoirTime2plot));


% For plotting of CO2 plumes
% bounds of 2008 plume:
ZoomX1 = 0.4375e6;
ZoomY1 = 6.47e6;
ZoomX2 = 0.4395e6;
ZoomY2 = 6.474e6;

figure; set(gcf, 'Position', [1 1 2000 600])
hold on

for i = 1:numel(ReservoirTime2plot)
    
    % get reservoir time index
    [rti,~] = find(sim_report.ReservoirTime==ReservoirTime2plot(i));

    % meaningful profiles
    densityCO2  = fluid.rhoG(states{rti}.pressure);  % fluid.rhoG is function handle to get CO2 density
    satCO2      = states{rti}.s(:,2);
    massCO2     = model.rock.poro.*model.G.cells.volumes.* model.G.cells.H.*satCO2.*densityCO2; % kg
    
    maxMassCO2(i)= max(massCO2);
    

    subplot(1,numel(ReservoirTime2plot),i)
    hold on

    % Add CO2 mass data: Note: To ensure cell data is plotted vertically
    % above the traps plotted by mapPlot, we modify the z-coordinate of the
    % faces to be z = -100
    Gt_tmp = Gt;
    Gt_tmp.nodes.z = -100*ones(Gt_tmp.nodes.num,1);
    plotFaces(Gt_tmp, boundaryFaces(Gt_tmp), 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt_tmp, massCO2, satCO2>CO2plumeOutline_SatTol, 'EdgeColor','none') % only plot plume that has sat > tolerance specified 
    title(num2str(Years2plot(i)), 'fontSize', 18); axis equal
    
    
    % Add all CO2 plume outlines that have a matching year to
    % Years2plot(i): Note: the plume outlines are plotted at with a
    % z-coordinate of -150 to ensure the outlines are on top of the plots
    % (for coloring purposes).
    for j = 1:numel(plume)
        if plume{j}.year == Years2plot(i)
            disp('Plotting Observed CO2 plume outline...')
            line(plume{j}.outline(:,1), plume{j}.outline(:,2), -150*ones(numel(plume{j}.outline(:,2)),1), 'LineWidth',3, 'Color','r')
        end
    end
    
    
    % Add injection point:
    % The following could be placed outside the subplot loop, or could be
    % plotted using mapPlot. Note: plot3() is used to ensure point is
    % plotted above other plots, with a z-coordinate of -200, so it remains
    % visible.
    
    % actual location
    plot3(wellXcoord, wellYcoord, -200, 'o', ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',10)
    % simulated location
    plot3(wellCoord_x, wellCoord_y, -200, 'x', ...
        'LineWidth',3,  ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',10)
    
    axis tight off
    box
    
    % We visualize the spill paths between structural traps.
    % Note: the plot is produced at z = 0.
    mapPlot(gcf, Gt, 'traps', trapstruct.traps, 'trapalpha', 0.2, ...
        'rivers', trapstruct.cell_lines, 'rivercolor', [1 0 0], ...
        'maplines', 20); % 'plumes',massCO2/1e9);
    
    % caution: adding this map changes values of colorbar, thus, it is
    % important to adjust the colorbar afterwards.
    
    % We then re-plot the contours of the grid topology, using contour3()
    % to plot the contours at the elevation determined by function inside
    % drawContours3(). To ensure these final contours are on top of all
    % other plots, we set cells.z to be negative values.
    ax = get(gcf, 'currentaxes');
    drawContours3(ax, Gt_tmp, -Gt_tmp.cells.z, 20, 'color', 'k');

end

% First adjust x and y limits (must do before adjusting colorbar!)
if ZoomIntoPlume
    set(findobj(gcf,'type','axes'),'xlim',[ZoomX1 ZoomX2]);
    set(findobj(gcf,'type','axes'),'ylim',[ZoomY1 ZoomY2]);
end

% Then adjust colorbar map. Works for both R2014a/earlier and later
% releases. The color bar is set between 0 and the maximum CO2 saturation
% (or mass) value that has occurred over all the years plotted.
cmax = max(maxMassCO2);
set(findobj(gcf,'type','axes'),'clim',[0, cmax]);

% Add a colorbar beside last subplot, and adjust it's position
[ hcb ] = setColorbarHandle( gcf, 'LabelName', 'kg', 'fontSize', 18 );
set(hcb,'Position',[0.915 0.185 0.01 0.662])


hfig = gcf;
hax  = gca;

end

