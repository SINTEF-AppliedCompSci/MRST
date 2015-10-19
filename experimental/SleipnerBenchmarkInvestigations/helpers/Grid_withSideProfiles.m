function [ hfig, hax ] = Grid_withSideProfiles( Gt, wellXcoord, wellYcoord, wellCoord_x, wellCoord_y, trapstruct, varargin )
% plot top view of grid
% plots slices taken through specified grid cell index(s)

    opt.figname = '';
    
    opt.plotNorthwardSlices = false; % if false, only east-west slices are plotted  
    opt.sliceCellIndex = [];
  
    opt = merge_options(opt, varargin{:});
    

    % Set-up number of figure subplots:
    figure('name',opt.figname); set(gcf,'Position',[1 1 1060 582]);
    numSlices = numel(opt.sliceCellIndex);
    if numSlices == 1
        if opt.plotNorthwardSlices
            sp1 = subplot(2,3,[1 4]);
            sp2 = subplot(2,3,2);
            sp2v = subplot(2,3,3);
        else
            sp1 = subplot(2,2,[1 3]);
            sp2 = subplot(2,2,2);
        end
        
    elseif numSlices == 2
        if opt.plotNorthwardSlices
            sp1 = subplot(2,3,[1 4]);
            sp2 = subplot(2,3,2);
            sp3 = subplot(2,3,5);
            sp2v = subplot(2,3,3);
            sp3v = subplot(2,3,6);
        else
            sp1 = subplot(2,2,[1 3]);
            sp2 = subplot(2,2,2);
            sp3 = subplot(2,2,4);
            %splabel = {'sp2', 'sp3'};
        end

    elseif numSlices == 3
        if opt.plotNorthwardSlices
            sp1 = subplot(3,3,[1 4 7]);
            sp2 = subplot(3,3,2);
            sp2v = subplot(3,3,3);
            sp3 = subplot(3,3,5);
            sp3v = subplot(3,3,6);
            sp4 = subplot(3,3,8);
            sp4v = subplot(3,3,9);
        else
            sp1 = subplot(3,2,[1 3 5]);
            sp2 = subplot(3,2,2);
            sp3 = subplot(3,2,4);
            sp4 = subplot(3,2,6);
            %splabel = {'sp2', 'sp3', 'sp4'};
        end

    else
        error('Only implemented for 2 or 3 slices')
    end
    


    %% Plot top view:
    subplot(sp1)
    hold on
    % Add CO2 mass data: Note: To ensure cell data is plotted vertically
    % above the traps plotted by mapPlot, we modify the z-coordinate of the
    % faces to be z = -100
    Gt_tmp = Gt;
    Gt_tmp.nodes.z = -100*ones(Gt_tmp.nodes.num,1);
    plotFaces(Gt_tmp, boundaryFaces(Gt_tmp), 'EdgeColor','k', 'LineWidth',3);
    %plotCellData(Gt_tmp, massCO2, satCO2>CO2plumeOutline_SatTol, 'EdgeColor','none'); % only plot plume that has sat > tolerance specified 
    %maxMassCO2 = max(massCO2);
    
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

    

    %% Plot slices:

    sliceLabels = ['A'; 'B'; 'C']; % ensure sliceLabels is as long as numSlices
    for sln = 1:numSlices
        
        sliceCellIndex_curr = opt.sliceCellIndex(sln);
        
        % Cross-sectional slices through a point:
        % first get index of point:
%         [ii,jj] = ind2sub(Gt.cartDims, sliceCellIndex_curr);
%         disp(['Slice cell index of ',num2str(sliceCellIndex_curr),' corresponds to I=',num2str(ii),', J=',num2str(jj)])
%         horiinx2color = find(Gt.cells.ij(:,2) == jj); % to color a horizontal line thru point
%         %vertinx2color = find(Gt.cells.ij(:,1) == ii); % to color a vertical line thru point
%         horiLineCoords = [Gt.cells.centroids(horiinx2color,1), Gt.cells.centroids(horiinx2color,2)];
%         %vertLineCoords = [Gt.cells.centroids(vertinx2color,1), Gt.cells.centroids(vertinx2color,2)];
%         horiLineCoords = [min(horiLineCoords); max(horiLineCoords)];

        ijk = gridLogicalIndices(Gt.parent);

        ii = ijk{1}(sliceCellIndex_curr);
        jj = ijk{2}(sliceCellIndex_curr);
        kk = ijk{3}(sliceCellIndex_curr);
        disp(['Slice cell index of ',num2str(sliceCellIndex_curr),' corresponds to I=',num2str(ii),', J=',num2str(jj)])
        
        horiinx2color = find(Gt.cells.ij(:,2) == jj); % to color a horizontal line thru point
        vertinx2color = find(Gt.cells.ij(:,1) == ii); % to color a vertical line thru point
        
        horiLineCoords = [Gt.cells.centroids(horiinx2color,1), Gt.cells.centroids(horiinx2color,2)];
        vertLineCoords = [Gt.cells.centroids(vertinx2color,1), Gt.cells.centroids(vertinx2color,2)];
        
        horiLineCoords = [min(horiLineCoords); max(horiLineCoords)];
        vertLineCoords = [min(vertLineCoords); max(vertLineCoords)];
        
        if opt.plotNorthwardSlices
            
            if numSlices == 1
                sp2 = subplot(2,3,2);

            elseif numSlices == 2
                if sln == 1
                    sp2 = subplot(2,3,2);
                elseif sln == 2
                    sp3 = subplot(2,3,5);
                end

            elseif numSlices == 3
                if sln == 1
                    sp2 = subplot(3,3,2);
                elseif sln == 2
                    sp3 = subplot(3,3,5);
                elseif sln == 3
                    sp4 = subplot(3,3,8);
                end
            end
            
        else
            
            if numSlices == 1
                sp2 = subplot(2,2,2);

            elseif numSlices == 2
                if sln == 1
                    sp2 = subplot(2,2,2);
                elseif sln == 2
                    sp3 = subplot(2,2,4);
                end

            elseif numSlices == 3
                if sln == 1
                    sp2 = subplot(3,2,2);
                elseif sln == 2
                    sp3 = subplot(3,2,4);
                elseif sln == 3
                    sp4 = subplot(3,2,6);
                end
            end
        
        end
        
        % plot grid of east-west slice
        plotGrid(Gt.parent, ijk{2} == jj, 'FaceColor', 'none'); view([0 0])
        
        % plot heights (max, free)
        %plotPlume(Gt.parent, Gt, states{rti}.h_max,  ijk{2} == jj, 'FaceColor',getInventoryColors(3),'EdgeColor','w','EdgeAlpha',.1)
        %plotPlume(Gt.parent, Gt, states{rti}.h_free, ijk{2} == jj, 'FaceColor',getInventoryColors(6),'EdgeColor','w','EdgeAlpha',.1)

        %xlabel('West to East','FontSize',16);
        %xlabel(['Slice ' sliceLabels(sln) ':' sliceLabels(sln)],'FontSize',16);
        %zlabel('depth (meters)','FontSize',16)
        %set(gca,'DataAspect',[1 1 1/50])
         %xlim([xmin xmax]);
         %zlim([zmin zmax]);
        %axis off;
        box;
        title(['slice ' sliceLabels(sln) ':' sliceLabels(sln)]);



        % Add slice line and label to TOPVIEW subplot Note: the slice line
        % is plotted using a z-coordinate of -250 such that the line is
        % visible ontop of the other previously plotted items.

        subplot( sp1 )
        hold on
        
        line(horiLineCoords(:,1), horiLineCoords(:,2), -250*ones(numel(horiLineCoords(:,2)),1), 'LineWidth',2, 'Color',[0.5 0.5 0.5], 'LineStyle','-.');

        % Add label at each end of line
        %text(horiLineCoords(1,1), horiLineCoords(1,2),   sliceLabels(sln), 'HorizontalAlignment','right');
        %text(horiLineCoords(2,1), horiLineCoords(2,2),   sliceLabels(sln), 'HorizontalAlignment','left');
        text(horiLineCoords(1,1)-7000, horiLineCoords(1,2),   sliceLabels(sln), 'HorizontalAlignment','left', 'FontSize',14);
        text(horiLineCoords(2,1)+7000, horiLineCoords(2,2),   sliceLabels(sln), 'HorizontalAlignment','right','FontSize',14);
    
        
        if opt.plotNorthwardSlices
           
            % plot grid of south-north slice
            if numSlices == 1
                sp2v = subplot(2,3,3);

            elseif numSlices == 2
                if sln == 1
                    sp2v = subplot(2,3,3);
                elseif sln == 2
                    sp3v = subplot(2,3,6);
                end

            elseif numSlices == 3
                if sln == 1
                    sp2v = subplot(3,3,3);
                elseif sln == 2
                    sp3v = subplot(3,3,6);
                elseif sln == 3
                    sp4v = subplot(3,3,9);
                end
            end

            plotGrid(Gt.parent, ijk{1} == ii, 'FaceColor', 'none'); view([90 0])
            box;
            title(['slice ' sliceLabels(sln) '_v:' sliceLabels(sln) '_v']);
            
            % Add slice line and label to TOPVIEW subplot
            subplot( sp1 )
            hold on
            line(vertLineCoords(:,1), vertLineCoords(:,2), -250*ones(numel(vertLineCoords(:,2)),1), 'LineWidth',2, 'Color',[0.5 0.5 0.5], 'LineStyle','-.');

            % Add label at each end of line
            text(vertLineCoords(1,1), vertLineCoords(1,2)-7000,   [sliceLabels(sln) '_v'], 'VerticalAlignment','bottom', 'FontSize',14);
            text(vertLineCoords(2,1), vertLineCoords(2,2)+7000,   [sliceLabels(sln) '_v'], 'VerticalAlignment','top',    'FontSize',14);

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

