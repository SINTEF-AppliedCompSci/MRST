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

    plotFaces(Gt, bf, 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt, massCO2/1e9, satCO2>CO2plumeOutline_SatTol, 'EdgeColor','none') % only plot plume that has sat > tolerance specified 
    %title({'Mass of CO2 at';['year ', num2str(Years2plot(i))]}, 'fontSize', 18); axis equal
    title(['year ', num2str(Years2plot(i))], 'fontSize', 18); axis equal
    
    % add colorbar if last subplot is being plotted: TODO --- make it look
    % better
    %if i == numel(ReservoirTime2plot)
        %%hcb = colorbar; hcb.Label.String = 'Mt'; set(hcb, 'fontSize', 18)
        [ ~ ] = setColorbarHandle( gcf, 'LabelName', 'Mt', 'fontSize', 18 );
    %end
    
    % add all CO2 plume outlines that have a matching year to Years2plot(i):
    for j = 1:numel(plume)
        if plume{j}.year == Years2plot(i)
            disp('Plotting Observed CO2 plume outline...')
            line(plume{j}.outline(:,1), plume{j}.outline(:,2), 'LineWidth',3, 'Color','r')
        end
    end
    
    % (The following could be placed outside the subplot loop)
    % Add injection point:
    % actual location
    plot(wellXcoord,wellYcoord,'o', ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',10)
    % simulated location
    plot(wellCoord_x,wellCoord_y,'x', ...
        'LineWidth',3,  ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',10)
    
    axis tight
    box
    
    % We visualize the spill paths between structural traps
    mapPlot(gcf, Gt, 'traps', trapstruct.traps, 'rivers', trapstruct.cell_lines, ...
        'maplines',20); % default maplines is 40
    % caution: adding this map changes values of colorbar (possibly to
    % elevation?). Thus, it is important to adjust the colorbar afterwards.

end

% First adjust x and y limits (must do before adjusting colorbar!)
if ZoomIntoPlume
    set(findobj(gcf,'type','axes'),'xlim',[ZoomX1 ZoomX2]);
    set(findobj(gcf,'type','axes'),'ylim',[ZoomY1 ZoomY2]);
end

% Then adjust colorbar map. Works for both R2014a/earlier and later releases.
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

