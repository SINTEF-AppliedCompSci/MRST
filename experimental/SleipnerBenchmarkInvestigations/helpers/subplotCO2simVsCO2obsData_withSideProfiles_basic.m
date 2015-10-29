function [ hfig, hax ] = subplotCO2simVsCO2obsData_withSideProfiles_basic( Years2plot, SimStartYear, plume, sim_report, Gt, states, fluid, model, wellXcoord, wellYcoord, wellCoord_x, wellCoord_y, trapstruct, CO2plumeOutline_SatTol, varargin )
% plot top view of grid with CO2 mass at a given states result
% plots slices taken through specified grid cell index(s)

    opt.figname = '';
    
    opt.ZoomIntoPlume = true;
    
    opt.sliceCellIndex_ew = []; % east-west slices
    opt.sliceCellIndex_sn = []; % south-north slices
    opt.SleipnerBounded = false; % if true, subplots will be bounded to set region.
    opt.legendWithFreeCO2Only = false;  % if true, what to include in legend
                                        % will be determined by comparing
                                        % free to residual heights.

    
    opt = merge_options(opt, varargin{:});
    
    % For plotting of CO2 plumes
    % bounds of 2008 plume:
    ZoomX1 = 0.4375e6;
    ZoomY1 = 6.47e6;
    ZoomX2 = 0.4395e6;
    ZoomY2 = 6.474e6;

    % Get limits of side view subplot:
    if opt.SleipnerBounded
        if opt.ZoomIntoPlume
            xmin = ZoomX1;
            xmax = ZoomX2;
            ymin = ZoomY1;
            ymax = ZoomY2;
        else
            xmin = 436914;
            xmax = 440114;
            ymin = 6469150;
            ymax = 6475050;
        end
        zmin = 790; %802.0627;
        zmax = 845; %850.6440;
    else
        xmin = min(Gt.nodes.coords(:,1));
        ymin = min(Gt.nodes.coords(:,2));
        xmax = max(Gt.nodes.coords(:,1));
        ymax = max(Gt.nodes.coords(:,2));
        zmax = max(Gt.cells.z+Gt.cells.H);       % the deepest depth
        zmin = min(Gt.cells.z);                  % the shallowest depth
    end

    % set up plot/subplot handles:
    % E-W slices will appear on same figure as Top View
    hfigEW = figure('name',opt.figname); set(gcf,'Position',[1 1 1060 582]);
    numSlices_ew = numel(opt.sliceCellIndex_ew);
    if numSlices_ew == 2
        % for 2 slices
        sp1 = subplot(2,2,[1 3]);
        sp2 = subplot(2,2,2);
        sp3 = subplot(2,2,4);

    elseif numSlices_ew == 3
        % for 3 slices
        sp1 = subplot(3,2,[1 3 5]);
        sp2 = subplot(3,2,2);
        sp3 = subplot(3,2,4);
        sp4 = subplot(3,2,6);

    else
        error('Only implemented for 2 or 3 east-west slices')
    end
    
    % S-N slices will appear on separate figure
    hfigSN = figure('name',[opt.figname 'SN']); set(gcf,'Position',[1 1 1060 582]);
    numSlices_sn = numel(opt.sliceCellIndex_sn);
    if numSlices_sn == 1
        % for 1 slice
        sp1_sn = subplot(2,2,[1 2]);

    elseif numSlices_sn == 2
        % for 2 slices
        sp1_sn = subplot(2,2,[1 2]);
        sp2_sn = subplot(2,2,[3 4]);

    else
        error('Only implemented for 1 or 2 south-north slices')
    end
    
    
    % get reservoir time index
    ReservoirTime2plot  = convertFrom( Years2plot - SimStartYear , year);   % seconds
    [rti,~]             = find(sim_report.ReservoirTime==ReservoirTime2plot);

    % meaningful profiles
    densityCO2  = fluid.rhoG(states{rti}.pressure);  % fluid.rhoG is function handle to get CO2 density
    satCO2      = states{rti}.s(:,2);
    massCO2     = model.rock.poro.*model.G.cells.volumes.* model.G.cells.H.*satCO2.*densityCO2; % kg
    
    % get the free plume height, and the max plume height reached in past
    if ~isempty(states)
        states = addCO2HeightData(states, Gt, fluid);
    end


    %% Plot top view:
    figure( hfigEW )
    subplot( sp1 )
    hold on
    % Add CO2 mass data: Note: To ensure cell data is plotted vertically
    % above the traps plotted by mapPlot, we modify the z-coordinate of the
    % faces to be z = -100
    Gt_tmp = Gt;
    Gt_tmp.nodes.z = -100*ones(Gt_tmp.nodes.num,1);
    plotFaces(Gt_tmp, boundaryFaces(Gt_tmp), 'EdgeColor','k', 'LineWidth',3);
    plotCellData(Gt_tmp, massCO2, satCO2>CO2plumeOutline_SatTol, 'EdgeColor','none'); % only plot plume that has sat > tolerance specified 
    maxMassCO2 = max(massCO2);
    
    % Add all CO2 plume outlines that have a matching year to
    % Years2plot(i): Note: the plume outlines are plotted at with a
    % z-coordinate of -150 to ensure the outlines are on top of the plots
    % (for coloring purposes).
    for j = 1:numel(plume)
        if plume{j}.year + 0.5 == Years2plot
            disp('Plotting Observed CO2 plume outline...')
            line(plume{j}.outline(:,1), plume{j}.outline(:,2), -150*ones(numel(plume{j}.outline(:,2)),1), 'LineWidth',3, 'Color','r')
            % adjust title: we assume plume data taken mid-year (i.e.,
            % July) thus we plot states results at 1999.5 and plume outline
            % at 1999.5, however title is left as 1999
            title(num2str(plume{j}.year), 'fontSize', 18);
        end
    end

    % Add injection point: The following could be placed outside the
    % subplot loop, or could be plotted using mapPlot. Note: plot3() is
    % used to ensure point is plotted above other plots, with a
    % z-coordinate of -200, so it remains visible.
    % actual location:
    plot3(wellXcoord, wellYcoord, -200, 'o', ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',10)
    % simulated location:
    plot3(wellCoord_x, wellCoord_y, -200, 'x', ...
        'LineWidth',3,  ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',10)

    % Adjust axis
    axis equal tight off


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


    % Adjustments to TOP VIEW subplot:
    % 1. adjust x and y limits (must do before adjusting colorbar!)
    if opt.ZoomIntoPlume
        set(findobj(gcf,'type','axes'),'xlim',[ZoomX1 ZoomX2]);
        set(findobj(gcf,'type','axes'),'ylim',[ZoomY1 ZoomY2]);
    end

    % 2. adjust colorbar map. Works for both R2014a/earlier and later
    % releases. The color bar is set between 0 and the maximum CO2 saturation
    % (or mass) value that has occurred over all the years plotted.
    cmax = max(maxMassCO2);
    set(findobj(gcf,'type','axes'),'clim',[0, cmax]);

    % 3. Add a colorbar beside last subplot, and adjust it's position
    [ hcb ] = setColorbarHandle( gcf, 'LabelName', 'kg', 'fontSize', 18 );


    %% Plot East-West slices:

    sliceLabels = ['A'; 'B'; 'C']; % ensure sliceLabels is as long as numSlices
    for sln = 1:numSlices_ew
        
        sliceCellIndex_curr = opt.sliceCellIndex_ew(sln);
        
        % Cross-sectional slices through a point:
        % first get index of point:
        [ii,jj] = ind2sub(Gt.cartDims, sliceCellIndex_curr);
        disp(['Slice cell index of ',num2str(sliceCellIndex_curr),' corresponds to I=',num2str(ii),', J=',num2str(jj)])
        
        % to color a horizontal line thru point
        horiinx2color = find(Gt.cells.ij(:,2) == jj); 
        horiLineCoords = [Gt.cells.centroids(horiinx2color,1), Gt.cells.centroids(horiinx2color,2)];
        horiLineCoords = [min(horiLineCoords); max(horiLineCoords)];

        ijk = gridLogicalIndices(Gt.parent);

        ii = ijk{1}(sliceCellIndex_curr);
        jj = ijk{2}(sliceCellIndex_curr);
        kk = ijk{3}(sliceCellIndex_curr);
        disp(['Slice cell index of ',num2str(sliceCellIndex_curr),' corresponds to I=',num2str(ii),', J=',num2str(jj)])

        if numSlices_ew == 2
            if sln == 1
                subplot( sp2 );
            elseif sln == 2
                subplot( sp3 );
            end
        elseif numSlices_ew == 3
            if sln == 1
                subplot( sp2 );
            elseif sln == 2
                subplot( sp3 );
            elseif sln == 3
                subplot( sp4 );
            end
        end
        
        % plot grid
        plotGrid(Gt.parent, ijk{2} == jj, 'FaceColor', 'none'); view([0 0])
        
        % plot heights (max, free)
        plotPlume(Gt.parent, Gt, states{rti}.h_max,  ijk{2} == jj, 'FaceColor',getInventoryColors(3),'EdgeColor','w','EdgeAlpha',.1)
        plotPlume(Gt.parent, Gt, states{rti}.h_free, ijk{2} == jj, 'FaceColor',getInventoryColors(6),'EdgeColor','w','EdgeAlpha',.1)

        xlim([xmin xmax]);
        zlim([zmin zmax]);
        box;
        title(['slice ' sliceLabels(sln) ':' sliceLabels(sln)]);

        % Add slice line and label to TOPVIEW subplot Note: the slice line
        % is plotted using a z-coordinate of -250 such that the line is
        % visible ontop of the other previously plotted items.
        subplot( sp1 )
        hold on
        
        line(horiLineCoords(:,1), horiLineCoords(:,2), -250*ones(numel(horiLineCoords(:,2)),1), 'LineWidth',2, 'Color',[0.5 0.5 0.5], 'LineStyle','-.');

        % Add label at each end of line
        text(xmin-200, horiLineCoords(1,2),   sliceLabels(sln,:), 'HorizontalAlignment','left', 'FontSize',14);
        text(xmax+200, horiLineCoords(2,2),   sliceLabels(sln,:), 'HorizontalAlignment','right','FontSize',14);
    
    end
    
    
    %% Plot south-north slices:
    
    %sliceLabels = ['AA'; 'BB'; 'CC']; % ensure sliceLabels is as long as numSlices
    sliceLabels = ['C'; 'D'];
    for sln = 1:numSlices_sn
        
        sliceCellIndex_curr = opt.sliceCellIndex_sn(sln);
        
        % Cross-sectional slices through a point:
        % first get index of point:
        [ii,jj] = ind2sub(Gt.cartDims, sliceCellIndex_curr);
        disp(['Slice cell index of ',num2str(sliceCellIndex_curr),' corresponds to I=',num2str(ii),', J=',num2str(jj)])
        
        % to color a vertical line thru point
        vertinx2color = find(Gt.cells.ij(:,1) == ii);
        vertLineCoords = [Gt.cells.centroids(vertinx2color,1), Gt.cells.centroids(vertinx2color,2)];
        vertLineCoords = [min(vertLineCoords); max(vertLineCoords)];

        ijk = gridLogicalIndices(Gt.parent);

        ii = ijk{1}(sliceCellIndex_curr);
        jj = ijk{2}(sliceCellIndex_curr);
        kk = ijk{3}(sliceCellIndex_curr);
        disp(['Slice cell index of ',num2str(sliceCellIndex_curr),' corresponds to I=',num2str(ii),', J=',num2str(jj)])

        figure( hfigSN );
        if numSlices_sn == 1
            if sln == 1
                subplot( sp1_sn );
            end
        elseif numSlices_ew == 2
            if sln == 1
                subplot( sp1_sn );
            elseif sln == 2
                subplot( sp2_sn );
            end
        end
        
        % plot grid
        plotGrid(Gt.parent, ijk{1} == ii, 'FaceColor', 'none'); view([90 0])
        
        % plot heights (max, free)
        plotPlume(Gt.parent, Gt, states{rti}.h_max,  ijk{1} == ii, 'FaceColor',getInventoryColors(3),'EdgeColor','w','EdgeAlpha',.1)
        plotPlume(Gt.parent, Gt, states{rti}.h_free, ijk{1} == ii, 'FaceColor',getInventoryColors(6),'EdgeColor','w','EdgeAlpha',.1)

        ylim([ymin ymax]);
        zlim([zmin zmax]);
        box;
        title(['slice ' sliceLabels(sln,:) ':' sliceLabels(sln,:)]);

        % Add slice line and label to TOPVIEW subplot Note: the slice line
        % is plotted using a z-coordinate of -250 such that the line is
        % visible ontop of the other previously plotted items.
        figure( hfigEW );
        subplot( sp1 )
        hold on
        
        line(vertLineCoords(:,1), vertLineCoords(:,2), -250*ones(numel(vertLineCoords(:,2)),1), 'LineWidth',2, 'Color',[0.5 0.5 0.5], 'LineStyle','-.');

        % Add label at each end of line
        text(vertLineCoords(1,1), ymin-250,  sliceLabels(sln,:), 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize',14);
        text(vertLineCoords(2,1), ymax+250,  sliceLabels(sln,:), 'HorizontalAlignment','center', 'VerticalAlignment','top',    'FontSize',14);
    
    end
    
    %% Final adjustments to plots
    
    figure( hfigSN )
    
    % Adjust fontsize of A---A, B---B, etc. font
    set(findobj(gcf,'type','Text'),'FontSize',18);
    
    % Adjust fontsize of slice titles
    set(findobj(gcf,'type','Axes'),'FontSize',14);
    
    
    figure( hfigEW )
    
    % Adjust fontsize of A---A, B---B, etc. font
    set(findobj(gcf,'type','Text'),'FontSize',18);
    
    % Adjust fontsize of slice titles
    set(findobj(gcf,'type','Axes'),'FontSize',14);

    
    %% Return variables:
    hfig = gcf;
    hax  = gca;



end

