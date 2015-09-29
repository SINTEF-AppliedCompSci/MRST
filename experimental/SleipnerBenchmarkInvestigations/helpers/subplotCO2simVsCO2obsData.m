function [ hfig, hax ] = subplotCO2simVsCO2obsData(Years2plot, inj_year, plume, sim_report, Gt, states, fluid, model, wellXcoord, wellYcoord, wellCoord_x, wellCoord_y, trapstruct, ZoomIntoPlume, ZoomX1, ZoomX2, ZoomY1, ZoomY2, CO2plumeOutline_SatTol)

% Reservoir Time2plot is an array of all times in seconds you wish to
% visualize the CO2 plume (saturation). If any times coorespond to the
% observation plume data, that plume outline is plotted as well.

ReservoirTime2plot  = (Years2plot - inj_year(1)+1 ).*(365*24*60*60); % seconds

bf = boundaryFaces(Gt);
maxMassCO2 = zeros(1,numel(ReservoirTime2plot));

figure; set(gcf, 'Position', [1 1 1500 1100])
hold on

for i = 1:numel(ReservoirTime2plot)
    
    % get reservoir time index
    [rti,~] = find(sim_report.ReservoirTime==ReservoirTime2plot(i));

    % meaningful profiles
    %press       = states{rti}.pressure;
    %pressDiffFromHydrostatic = press - initState.pressure;
    densityCO2  = fluid.rhoG(states{rti}.pressure);  % fluid.rhoG is function handle to get CO2 density
    satCO2      = states{rti}.s(:,2);
    massCO2     = model.rock.poro.*model.G.cells.volumes.* model.G.cells.H.*satCO2.*densityCO2; % kg
    
    maxMassCO2(i)= max(massCO2);
    
    if numel(ReservoirTime2plot)==6
        subplot(2,3,i)
    elseif numel(ReservoirTime2plot)==3
        subplot(1,3,i)
    end
    %subplot(2, numel(ReservoirTime2plot)/2, i)
    hold on

    % Add CO2 mass data: Note: To ensure cell data is plotted vertically
    % above the traps plotted by mapPlot, we modify the z-coordinate of the
    % faces to be z = -100
    Gt_tmp = Gt;
    Gt_tmp.nodes.z = -100*ones(Gt_tmp.nodes.num,1);
    plotFaces(Gt_tmp, boundaryFaces(Gt_tmp), 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt_tmp, massCO2/1e9, satCO2>CO2plumeOutline_SatTol, 'EdgeColor','none') % only plot plume that has sat > tolerance specified 
    %title({'Mass of CO2 at';['year ', num2str(Years2plot(i))]}, 'fontSize', 18); axis equal
    title(['year ', num2str(Years2plot(i))], 'fontSize', 18); axis equal
    
    
    % Add colorbar (TODO: if last subplot is being plotted)
    %if i == numel(ReservoirTime2plot)
        %%hcb = colorbar; hcb.Label.String = 'Mt'; set(hcb, 'fontSize', 18)
        [ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Mt', 'fontSize', 18 );
    %end
    
    
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
    
    axis tight
    box
    
    % We visualize the spill paths between structural traps.
    % Note: the plot is produced at z = 0.
    mapPlot(gcf, Gt, 'traps', trapstruct.traps, 'rivers', trapstruct.cell_lines, ...
        'maplines',20, 'trapalpha',0.2); % 'plumes',massCO2/1e9);
    
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
cmax = max(maxMassCO2./1e9);
set(findobj(gcf,'type','axes'),'clim',[0, cmax]);


% % only works for R2014a/earlier
% hfig = gcf;
% subplots = get(hfig,'Children');
% 
% cmax = max(maxMassCO2./1e9);
% for i=1:length(subplots)
%     caxis(subplots(i),[0,cmax]);
% end
% 
% % only works for later than R2014a
% % make plotting adjustments to subplots
% axesHandles = get(gcf,'children');
% 
% % set caxis to be between 0 and the max CO2 mass value plotted in any of
% % the subplots
% cmax = max(maxMassCO2./1e9);
% for i=1:numel(axesHandles)
%     if strcmpi(axesHandles(i).Type,'axes')
%         axesHandles(i).CLim = [0 cmax];
%         
%         if ZoomIntoPlume
%            axesHandles(i).XLim = [ZoomX1 ZoomX2];
%            axesHandles(i).YLim = [ZoomY1 ZoomY2];
%         end
%     end
% end


hfig = gcf;
hax  = gca;

end

