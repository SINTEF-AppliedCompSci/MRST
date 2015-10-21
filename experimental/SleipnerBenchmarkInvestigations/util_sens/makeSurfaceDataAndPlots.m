function [plumes, topsurface, topfit, hCO2] = makeSurfaceDataAndPlots(plumes, Gt, varargin)

% Gt        - top surface grid
% plumes    - XY coordinates of plume outline

% interpolation is done such that plume and Grid top surface have
% intersecting XY coordinates:
%   topsurface    - returns z coordinate at (x,y) of the grid's caprock
%   topfit        - returns z coordinate at (x,y) of a planar surface
%                   (planar surface is a fit to the plume outline and to
%                   the grid's caprock).

opt.plotsOn = false;

opt = merge_options(opt, varargin{:});


%% Function handle of top-surface
if(isfield(Gt,'cartDims') && Gt.cells.num==prod(Gt.cartDims))
    X = reshape(Gt.cells.centroids(:,1),Gt.cartDims(1),Gt.cartDims(2));
    Y = reshape(Gt.cells.centroids(:,2),Gt.cartDims(1),Gt.cartDims(2));
    Z = reshape(Gt.cells.z,Gt.cartDims(1),Gt.cartDims(2));
    var.topsurface = @(coord) interp2(X',Y',Z',coord(:,1),coord(:,2));
else
    F = scatteredInterpolant(Gt.cells.centroids(:,1),Gt.cells.centroids(:,2),Gt.cells.z);
    var.topsurface = @(coord) F(coord(:,1),coord(:,2));
end

%% Function handle of top-fit (unique fit for each plume polygon)
for i=1:numel(plumes)
    line_coord = plumes{i}.outline;
    
    nd  = size(line_coord,1);
    A   = [ones(nd,1),line_coord];
    rhs = var.topsurface(line_coord);
    vec = A\rhs;
    var.topfit = @(coord) [ones(size(coord,1),1),coord]*vec;
    
    hCO2_tmp    = @(coord)(var.topfit(coord)-var.topsurface(coord));
    hCO2        = @(coord) hCO2_tmp(coord).*(hCO2_tmp(coord)>0).*insidePolygon(line_coord,coord);
    plumes{i}.hCO2  = hCO2;
    plumes{i}.h     = hCO2(Gt.cells.centroids);
     
end

%% Generate some plots (optional)
if opt.plotsOn
    [insidePlumeXY_tmp, insidePlumeXY, insidePlumeXYZ] = getInsidePlumeCoords(Gt, line_coord);
    [ ~ ] = makePlot1(Gt, line_coord);
    [ ~ ] = makePlot2(Gt, line_coord);    
    [ ~ ] = makePlot3(Gt, line_coord); 
    [ ~ ] = makePlot4(Gt, line_coord, 'SleipnerBounded',true, 'plume2008Bounded',false);
    [ ~ ] = makePlot5(Gt, hCO2, line_coord, insidePlumeXY);
    [ ~ ] = makePlot6(Gt, hCO2, line_coord);
end

%% Return variables
topsurface = var.topsurface;
topfit = var.topfit;
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%     Helpers functions - for plotting   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [insidePlumeXY_tmp, insidePlumeXY, insidePlumeXYZ] = getInsidePlumeCoords(Gt, line_coord)

    % then get XY coordinates of plume inside
    insidePlumeXY_tmp = [Gt.cells.centroids(:,1).*insidePolygon(line_coord, Gt.cells.centroids), ...
                         Gt.cells.centroids(:,2).*insidePolygon(line_coord, Gt.cells.centroids)];
    tmp1 = insidePlumeXY_tmp(:,1);
    tmp2 = insidePlumeXY_tmp(:,2);
    insidePlumeXY = [tmp1(insidePlumeXY_tmp(:,1)~=0), tmp2(insidePlumeXY_tmp(:,2)~=0)];

    % XYZ coordinates of plume inside
    insidePlumeXYZ = [insidePlumeXY_tmp, var.topsurface(insidePlumeXY_tmp)];
end

% -------------------------------------------------------------------------

function [ hfig ] = makePlot1(Gt, line_coord)

    figure; set(gcf,'Position',[1 1 800 700])

    % first plot top surface of grid
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'k', 'EdgeAlpha', .1);

    % adjust plot
    axis equal tight; view(35,25)
    set(gca,'DataAspect',[1 1 0.02]);
    xlabel('x','Fontweight','bold');
    ylabel('y','Fontweight','bold');
    zlabel('depth (meters)','Fontweight','bold');
    box;
    hold on


    % plot XYZ plume coordinates, using Z from top surface of grid
    offSet = 0; % meters
    hp1 = plot3(line_coord(:,1), line_coord(:,2), var.topsurface(line_coord)-offSet, 'r', 'LineWidth',2);

    % then color cells of Gt that are inside plume
    %plotCellData(Gt, insidePlumeXYZ, insidePlumeXYZ(:,1)~=0, 'FaceColor','r')

    % plot XYZ plume coordinates, using Z from the fitted planar surface
    hp2 = plot3(line_coord(:,1), line_coord(:,2), var.topfit(line_coord)-offSet, 'b', 'LineWidth',2);

    % adjust axis fontsize
    set(gca,'FontSize',14)

    % add legend
    hl = legend([hp1, hp2],{'plume outline following grid topography','plume outline with planar shoreline'},'Location','ne');
    set(hl,'FontSize',14, 'Location','northwest')

    hfig = gcf;

end

% -------------------------------------------------------------------------

function [ hfig ] = makePlot2(Gt, line_coord)

    % simple side view (facing west) without grid
    
    figure; set(gcf,'Position',[1 1 650 450])
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'none');%, 'EdgeAlpha', .1);
    hold on
    %plot(line_coord(:,2), topsurface(line_coord), 'r', 'LineWidth',2);
    %plot(line_coord(:,2), topfit(line_coord), 'b', 'LineWidth',2);

    % plot XYZ plume coordinates, using Z from top surface of grid
    %plot3(insidePlumeXY(:,1), insidePlumeXY(:,2), topsurface(insidePlumeXY), 'xr')
    offSet = 0; % meters
    hp1 = plot3(line_coord(:,1), line_coord(:,2), var.topsurface(line_coord)-offSet, 'r', 'LineWidth',2);

    % then color cells of Gt that are inside plume
    %plotCellData(Gt, insidePlumeXYZ, insidePlumeXYZ(:,1)~=0, 'FaceColor','r')

    % plot XYZ plume coordinates, using Z from the fitted planar surface
    %plot3(insidePlumeXY(:,1), insidePlumeXY(:,2), topfit(insidePlumeXY), 'sqb')
    hp2 = plot3(line_coord(:,1), line_coord(:,2), var.topfit(line_coord)-offSet, 'b', 'LineWidth',2);

    view([90,0])
    ylim([6470000 6474500])
    zlim([800 820])

    % adjust axis fontsize
    set(gca,'FontSize',14)
    zlabel('depth (meters)')
    ylabel('South to North')

    % add legend
    hl = legend([hp1, hp2],{'plume outline following grid topography','plume outline with planar shoreline'},'Location','nw');
    set(hl,'FontSize',16)

    box
    hfig = gcf; 
end

% -------------------------------------------------------------------------

function [ hfig ] = makePlot3(Gt, line_coord)

    % simple side view (facing west)
    figure; set(gcf,'Position',[1 1 1200 400])

    % first plot top surface of grid
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'k', 'EdgeAlpha', .1);

    % adjust plot
    axis equal tight; view(35,25)
    set(gca,'DataAspect',[1 1 0.02]);
    xlabel('x','FontSize',16,'FontWeight','bold');
    ylabel('South to North','FontSize',16,'FontWeight','bold');
    zlabel('depth (meters)','FontSize',16,'FontWeight','bold');
    box;

    hold on


    % plot XYZ plume coordinates, using Z from top surface of grid
    %plot3(insidePlumeXY(:,1), insidePlumeXY(:,2), topsurface(insidePlumeXY), 'xr')
    offSet = 0; % meters
    hp1 = plot3(line_coord(:,1), line_coord(:,2), var.topsurface(line_coord)-offSet, 'r', 'LineWidth',2);

    % then color cells of Gt that are inside plume
    %plotCellData(Gt, insidePlumeXYZ, insidePlumeXYZ(:,1)~=0, 'FaceColor','r')

    % plot XYZ plume coordinates, using Z from the fitted planar surface
    %plot3(insidePlumeXY(:,1), insidePlumeXY(:,2), topfit(insidePlumeXY), 'sqb')
    hp2 = plot3(line_coord(:,1), line_coord(:,2), var.topfit(line_coord)-offSet, 'b', 'LineWidth',2);

    view([90,0])
    ylim([6470000 6474500])
    zlim([800 825])

    % adjust axis fontsize
    set(gca,'FontSize',14)

    % add legend
    hl = legend([hp1, hp2],{'plume outline following grid topography','plume outline with planar shoreline'},'Location','n');
    set(hl,'FontSize',16)
    
    hfig = gcf;

end

% -------------------------------------------------------------------------

function [ hfig ] = makePlot4(Gt, line_coord, varargin)


    opt.SleipnerBounded = true;
    opt.plume2008Bounded = false;

    opt = merge_options(opt, varargin{:});

    % side view figures    
    figure; set(gcf,'Position',[1 1 1700 350]);

    % Set up subplot sizes:
    if opt.SleipnerBounded
        xmin = 436914;
        xmax = 440114;
        ymin = 6469150;
        ymax = 6475050;
        zmin = 790; %802.0627;
        zmax = 845; %850.6440;
    elseif opt.plume2008Bounded
        % bounds of 2008 plume:
        ZoomX1 = 0.438e6; %0.4375e6;
        ZoomY1 = 6.4705e6; %6.47e6;
        ZoomX2 = 0.4395e6;
        ZoomY2 = 6.474e6;
        xmin = ZoomX1; xmax = ZoomX2;
        ymin = ZoomY1; ymax = ZoomY2;
        zmax = 825;
        zmin = 800;
    else
        xmin = min(Gt.nodes.coords(:,1));
        ymin = min(Gt.nodes.coords(:,2));
        xmax = max(Gt.nodes.coords(:,1));
        ymax = max(Gt.nodes.coords(:,2));
        zmax = max(Gt.cells.z+Gt.cells.H);       % the deepest depth
        zmin = min(Gt.cells.z);                  % the shallowest depth
    end
    width1 = xmax - xmin;
    width2 = ymax - ymin;
    scaleFactor = width1/width2;
    normWidth1 = 0.3;
    normWidth2 = normWidth1/scaleFactor;

    left1 = 0.07;
    bottom = 0.17;
    gap = 0.06;
    normHeight = 0.75;
    hsp1 = subplot('Position', [left1                bottom normWidth1 normHeight]);
    hsp2 = subplot('Position', [left1+normWidth1+gap bottom normWidth2 normHeight]);


    subplot(hsp1)
    % grid
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'k', 'EdgeAlpha', .1);
    % plume outline on top of grid
    hold on;
    plot3(line_coord(:,1), line_coord(:,2), var.topsurface(line_coord), 'r', 'LineWidth',2)
    plot3(line_coord(:,1), line_coord(:,2), var.topfit(line_coord), 'b', 'LineWidth',2)

    view([0,0])
    xlabel('West to East','FontSize',16);
    zlabel('depth (meters)','FontSize',16)
    %set(gca,'DataAspect',[1 1 1/50])
    xlim([xmin xmax]);
    zlim([zmin zmax]);
    box
    set(gca,'FontSize',14)

    subplot(hsp2)
    % grid
    plotGrid(Gt, 'FaceColor', 'none', 'EdgeColor', 'k', 'EdgeAlpha', .1);
    % plume outline on top of grid
    hold on;
    plot3(line_coord(:,1), line_coord(:,2), var.topsurface(line_coord), 'r', 'LineWidth',2)
    plot3(line_coord(:,1), line_coord(:,2), var.topfit(line_coord), 'b', 'LineWidth',2)

    view([90,0])
    ylabel('South to North','FontSize',16);
    zlabel('depth (meters)','FontSize',16)
    %set(gca,'DataAspect',[1 1 1/50])
    ylim([ymin ymax]);
    zlim([zmin zmax]);
    box
    set(gca,'FontSize',14)

    hfig = gcf; 

end

% -------------------------------------------------------------------------

function [ hfig ] = makePlot5(Gt, hCO2, line_coord, insidePlumeXY)


    figure; set(gcf,'Position',[1 1 1300 500])
    subplot(1,2,1)

    hold on
    plot3(insidePlumeXY(:,1), insidePlumeXY(:,2), var.topsurface(insidePlumeXY), 'xr')
    plot3(insidePlumeXY(:,1), insidePlumeXY(:,2), var.topfit(insidePlumeXY), 'sqb')
    plot3(line_coord(:,1), line_coord(:,2), var.topsurface(line_coord), 'r', 'LineWidth',2)
    plot3(line_coord(:,1), line_coord(:,2), var.topfit(line_coord), 'b', 'LineWidth',2)
    view(3)

    %title({['Elevations of (x,y) inside plume as determined by'];['functions (wrt ',num2str(plumes{i}.year),' plume).']})
    legend('topsurface(x,y)','topfit(x,y)','topsurface(x_p,y_p)','topfit(x_p,y_p)','Location','NorthEast')
    axis equal tight;
    set(gca,'DataAspect',[1 1 1/50])
    grid
    zlabel('depth (meters)')

    subplot(1,2,2)
    plotCellData(Gt, hCO2(Gt.cells.centroids)); colorbar
    line(line_coord(:,1), line_coord(:,2), var.topsurface(line_coord)-1, 'LineWidth',2, 'Color','r'); % -1 to ensure visibility of plume outline against other plots
    view(2)

    title('hCO2: (topfit - topsurface) > 0, inside plume')
    legend('hCO2(x,y)','topsurface(x_p,y_p)', 'Location','NorthEast')
    axis equal tight;
    set(gca,'DataAspect',[1 1 1/50])
    grid
    zlabel('depth (meters)')
    
    hfig = gcf;

end

% -------------------------------------------------------------------------

function [ hfig ] = makePlot6(Gt, hCO2, line_coord)

    % show observed CO2 heights (that were computed using a grid and a
    % plume outline)
    figure; %set(gcf,'Position',[1 1 500 500])

    plotCellData(Gt, hCO2(Gt.cells.centroids), 'EdgeColor','none'); colorbar
    %line(line_coord(:,1), line_coord(:,2), var.topsurface(line_coord)-3, 'LineWidth',2, 'Color','r'); % -1 to ensure visibility of plume outline against other plots
    view(2)

    %title('hCO2: (topfit - topsurface) > 0, inside plume')
    title('CO_2 heights (between caprock and planar shoreline)')
    %legend('hCO2(x,y)','topsurface(x_p,y_p)', 'Location','NorthEast')
    axis equal tight;
    set(gca,'DataAspect',[1 1 1/50])
    grid
    zlabel('depth (meters)')
    
    hfig = gcf;

end


%%%%%%%%%%%%%%%%%%         End of Helper Functions        %%%%%%%%%%%%%%%%%


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%         Other Helper Functions         %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

