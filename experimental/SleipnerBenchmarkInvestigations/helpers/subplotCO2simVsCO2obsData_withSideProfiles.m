function [ hfig, hax ] = subplotCO2simVsCO2obsData_withSideProfiles( Years2plot, SimStartYear, plume, sim_report, Gt, states, fluid, model, wellXcoord, wellYcoord, wellCoord_x, wellCoord_y, trapstruct, CO2plumeOutline_SatTol, varargin )
% new plotting routine to include side profiles slices for each subplot


    opt.figname = '';
    
    opt.ZoomIntoPlume = true;
    
    opt.sliceCellIndex = [];
    opt.SleipnerBounded = false; % if true, subplots will be bounded to set region.
    opt.legendWithFreeCO2Only = false;  % if true, what to include in legend
                                        % will be determined by comparing
                                        % free to residual heights.

    
    opt = merge_options(opt, varargin{:});
    
    %% Set up figure:
    % get correct number of subplots, depending on the number of slice cell
    % index that was requested (if any).
    if isempty(opt.sliceCellIndex)
        opt.sliceCellIndex = model.wellmodel.W.cells;
    end
    numSlices       = numel(opt.sliceCellIndex);
    % to ensure a top view subplot gets plotted in at least twice as many subplots as a slice
    if numSlices == 1
        numRows4TopView = 2;
    elseif numSlices == 2
        numRows4TopView = 2;
    else
        numRows4TopView = ceil(numSlices/2);
    end
    numPlotRows     = numSlices + numRows4TopView;
    numPlotCols     = numel(Years2plot);
    addRow = 1; % number of rows to place inbetween top view and side profile rows

    ReservoirTime2plot  = convertFrom( Years2plot - SimStartYear , year); % seconds

    maxMassCO2 = zeros(1,numel(ReservoirTime2plot));


    % For plotting of CO2 plumes
    % bounds of 2008 plume:
    ZoomX1 = 0.4375e6;
    ZoomY1 = 6.47e6;
    ZoomX2 = 0.4395e6;
    ZoomY2 = 6.474e6;

    
    figure('name',opt.figname);
    if numPlotCols == 6
        set(gcf, 'Position', [1 1 2000 900])
    elseif numPlotCols == 1
        set(gcf, 'Position', [1 1 500 1000])
    else
        disp('You will need to adjust the figure window size.')
    end
    hold on

    %% First plot all subplots of CO2 mass with plume outlines from top view:
    for i = 1:numel(ReservoirTime2plot)

        % get reservoir time index
        [rti,~] = find(sim_report.ReservoirTime==ReservoirTime2plot(i));

        % meaningful profiles
        densityCO2  = fluid.rhoG(states{rti}.pressure);  % fluid.rhoG is function handle to get CO2 density
        satCO2      = states{rti}.s(:,2);
        massCO2     = model.rock.poro.*model.G.cells.volumes.* model.G.cells.H.*satCO2.*densityCO2; % kg

        maxMassCO2(i)= max(massCO2);


        % To plot Top Views in bottom 'N=numRows4TopViewrow' rows subplots:
        %subplot( numPlotRows, numPlotCols, [i+(numPlotRows-1)*numPlotCols i+numRows4TopView*numPlotCols] )
        
        % To plot Top Views in first 'N=numRows4TopViewrow' rows subplots:
        subplot( numPlotRows+addRow, numPlotCols, [i i+numRows4TopView*numPlotCols] )

        hold on

        % Add CO2 mass data: Note: To ensure cell data is plotted
        % vertically above the traps plotted by mapPlot, we modify the
        % z-coordinate of the faces to be z = -100
        Gt_tmp = Gt;
        Gt_tmp.nodes.z = -100*ones(Gt_tmp.nodes.num,1);
        plotFaces(Gt_tmp, boundaryFaces(Gt_tmp), 'EdgeColor','k', 'LineWidth',3);
        plotCellData(Gt_tmp, massCO2, satCO2>CO2plumeOutline_SatTol, 'EdgeColor','none') % only plot plume that has sat > tolerance specified 
        title(num2str(Years2plot(i)), 'fontSize', 18);
        axis equal


        % Add all CO2 plume outlines that have a matching year to
        % Years2plot(i): Note: the plume outlines are plotted at with a
        % z-coordinate of -150 to ensure the outlines are on top of the
        % plots (for coloring purposes).
        for j = 1:numel(plume)
            if plume{j}.year + 0.5 == Years2plot(i)
                disp('Plotting Observed CO2 plume outline...')
                line(plume{j}.outline(:,1), plume{j}.outline(:,2), -150*ones(numel(plume{j}.outline(:,2)),1), 'LineWidth',3, 'Color','r')
                % adjust title: we assume plume data taken mid-year (i.e.,
                % July) thus we plot states results at 1999.5 and plume
                % outline at 199.5, however title is left as 1999
                title(num2str(plume{j}.year), 'fontSize', 18);
            end
        end

        % Add injection point: The following could be placed outside the
        % subplot loop, or could be plotted using mapPlot. Note: plot3() is
        % used to ensure point is plotted above other plots, with a
        % z-coordinate of -200, so it remains visible.

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

        % We then re-plot the contours of the grid topology, using
        % contour3() to plot the contours at the elevation determined by
        % function inside drawContours3(). To ensure these final contours
        % are on top of all other plots, we set cells.z to be negative
        % values.
        ax = get(gcf, 'currentaxes');
        drawContours3(ax, Gt_tmp, -Gt_tmp.cells.z, 20, 'color', 'k');
        

    end

    % Adjustments to TOP VIEW subplots:
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
    if numPlotCols == 6
        set(hcb,'Position',[0.915 0.185 0.01 0.662])
    elseif numPlotCols == 1
        set(hcb,'Position',[0.78 0.535 0.028 0.392])
    else
        disp('You will need to adjust the color bar position.')
    end
    



    %% Second plot all subplots of CO2 heights side profiles:
    
    % if no slice cell index was specified, default is to slice through the
    % injection point. Otherwise, 1 or more slice cell index can be
    % specified. If more than 1 slice is requested, figure subplots must be
    % set up to handle the correct number of TOP VIEW and SLICE plots.
    
    
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
    
    
    % get the free plume height, and the max plume height reached in past
    if ~isempty(states)
        states = addCO2HeightData(states, Gt, fluid);
    end
    
    sliceLabels = ['A'; 'B'; 'C'; 'D'; 'E']; % ensure sliceLabels is as long as numSlices
    for sln = 1:numSlices
        
        sliceCellIndex_curr = opt.sliceCellIndex(sln);
        
        % Cross-sectional slices through a point:
        % first get index of point:
        [ii,jj] = ind2sub(Gt.cartDims, sliceCellIndex_curr);
        disp(['Slice cell index of ',num2str(sliceCellIndex_curr),' corresponds to I=',num2str(ii),', J=',num2str(jj)])
        horiinx2color = find(Gt.cells.ij(:,2) == jj); % to color a horizontal line thru point
        %vertinx2color = find(Gt.cells.ij(:,1) == ii); % to color a vertical line thru point
        horiLineCoords = [Gt.cells.centroids(horiinx2color,1), Gt.cells.centroids(horiinx2color,2)];
        %vertLineCoords = [Gt.cells.centroids(vertinx2color,1), Gt.cells.centroids(vertinx2color,2)];
        horiLineCoords = [min(horiLineCoords); max(horiLineCoords)];

        ijk = gridLogicalIndices(Gt.parent);

        ii = ijk{1}(sliceCellIndex_curr);
        jj = ijk{2}(sliceCellIndex_curr);
        kk = ijk{3}(sliceCellIndex_curr);
        disp(['Slice cell index of ',num2str(sliceCellIndex_curr),' corresponds to I=',num2str(ii),', J=',num2str(jj)])

        
        for i = 1:numel(ReservoirTime2plot)

            % get reservoir time index
            [rti,~] = find(sim_report.ReservoirTime==ReservoirTime2plot(i));

            
            % If side views are in first rows of subplots:
            %subplot(numPlotRows, numPlotCols, i + (sln-1)*numPlotCols)
            
            % If side views come after rows used for Top Views:
            subplot( numPlotRows+addRow, numPlotCols, (i + numRows4TopView*numPlotCols + (sln-1)*numPlotCols) + addRow*numPlotCols)
            
            hold on

            % plot grid
            plotGrid(Gt.parent, ijk{2} == jj, 'FaceColor', 'none'); view([0 0])

            % plot heights (max, free)
            plotPlume(Gt.parent, Gt, states{rti}.h_max,  ijk{2} == jj, 'FaceColor',getInventoryColors(3),'EdgeColor','w','EdgeAlpha',.1)
            plotPlume(Gt.parent, Gt, states{rti}.h_free, ijk{2} == jj, 'FaceColor',getInventoryColors(6),'EdgeColor','w','EdgeAlpha',.1)

            %xlabel('West to East','FontSize',16);
            %xlabel(['Slice ' sliceLabels(sln) ':' sliceLabels(sln)],'FontSize',16);
            %zlabel('depth (meters)','FontSize',16)
            %set(gca,'DataAspect',[1 1 1/50])
            xlim([xmin xmax]);
            zlim([zmin zmax]);
            axis off;
            title(['slice ' sliceLabels(sln) ':' sliceLabels(sln)]);



            % Add slice line and label to TOPVIEW subplot
            % Note: the slice line is plotted using a z-coordinate of -250 such
            % that the line is visible ontop of the other previously plotted
            % items.
            % If Top Views are in bottom 'N=numRows4TopViewrow' rows subplots:
            %subplot( numPlotRows, numPlotCols, [i+(numPlotRows-1)*numPlotCols i+numRows4TopView*numPlotCols] )
        
            % If Top Views are in first 'N=numRows4TopViewrow' rows subplots:
            subplot( numPlotRows+addRow, numPlotCols, [i i+numRows4TopView*numPlotCols] )
            
            hold on
            %plotCellData(Gt, -250*ones(Gt.cells.num,1), horiinx2color, 'FaceColor','red', 'EdgeColor','none')
            %plotGrid(Gt, 'FaceColor','none');
            line(horiLineCoords(:,1), horiLineCoords(:,2), -250*ones(numel(horiLineCoords(:,2)),1), 'LineWidth',2, 'Color',[0.5 0.5 0.5], 'LineStyle','-.');

            % Add label at each end of line
            %text(horiLineCoords(1,1), horiLineCoords(1,2),   sliceLabels(sln), 'HorizontalAlignment','right');
            %text(horiLineCoords(2,1), horiLineCoords(2,2),   sliceLabels(sln), 'HorizontalAlignment','left');
            text(xmin-200, horiLineCoords(1,2),   sliceLabels(sln), 'HorizontalAlignment','left', 'FontSize',14);
            text(xmax+200, horiLineCoords(2,2),   sliceLabels(sln), 'HorizontalAlignment','right','FontSize',14);
        
        end
    
    end
    
    % Adjust fontsize of A---A, B---B, etc. font
    set(findobj(gcf,'type','Text'),'FontSize',18);
    
    % Adjust fontsize of slice titles
    set(findobj(gcf,'type','Axes'),'FontSize',14);



    %% Return variables:
    hfig = gcf;
    hax  = gca;


end

