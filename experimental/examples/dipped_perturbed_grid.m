function Gt = dipped_perturbed_grid(varargin)
% Construct grid used in /co2lab/examples/papers/CAGEO-75/trappingExample1.m
%
% nx, ny - grid discretization (100x100 or 50x25 etc)
% Lx, Ly, H - dimensions of grid (in meters)
% Grid is 10 km by 5 km by 50 meters deep. 
%

    [opt.Lx, opt.Ly, opt.H] = deal(10000, 5000, 50);
    [opt.nx, opt.ny] = deal(100, 50);
    [opt.dz_s, opt.dz_n, opt.dz_e] = deal(0, 0, 0);
    opt.originalGrid = true;
    
    opt = merge_options(opt, varargin{:});
    
    G  = cartGrid([opt.nx opt.ny 1],[opt.Lx opt.Ly opt.H]);
    x  = G.nodes.coords(1:G.nodes.num/2,1)/opt.Lx; 
    y  = G.nodes.coords(1:G.nodes.num/2,2)/opt.Ly;
    z  = G.nodes.coords(1:G.nodes.num/2,3)/opt.H;
    if opt.originalGrid
        zt = z + x - 0.2*sin(5*pi*x).*sin(5*pi*y.^1.5) - 0.075*sin(1.25*pi*y) + 0.15*sin(x+y);
    else
        % ALPHA impacts amplitude of perturbed waves on top surface
        % WEDGE impacts aspect ratio of domain as a wedge (from north to south)
        % BEND impacts curve of top surface (from east to west)
        % SLOPE impacts top surface's dipping angle
        ALPHA = 0.2;
        WEDGE = 0.15; %0.15; 
        BEND = 0.075;
        SLOPE = 1;
        gamma1 = 5;
        gamma2 = 5;
        beta = 1.5; %1.5; % impacts amount of squeeze/stretch in waves
        phi = 1.25; %1.25; % impacts curvature of top from north to south

        fac = 2;
        xx  = min(max(0,(x-fac/opt.nx)/(1-fac*2/opt.nx)),1);
        %xx = x;
        yy  = min(max(0,(y-fac/opt.ny)/(1-fac*2/opt.ny)),1);
        %yy = y;
        oss = ALPHA*sin(gamma1*pi*xx).*sin(gamma2*pi*yy.^beta);
        %pad = gamma1*pi*x<pi*fac | gamma1*pi*(1-x)<pi*fac | gamma2*pi*y<pi*fac | gamma2*pi*(1-y)<pi*fac; 
        zt = z + SLOPE*x - oss - BEND*sin(phi*pi*y) +  WEDGE*sin(x+y);%-oss.*(~pad);
    end
    zb = 1 + x;
    %zt = zb-1;
    G.nodes.coords(:,3) = [zt; zb]*opt.H+1000;
    G = computeGeometry(G);
    
    
    % Alter the external boundaries to ensure trap regions do not touch
    % edges (i.e., forces co2 to spill out along these sides)
    finx_sb = boundaryFaceIndices(G,'South');
    finx_nb = boundaryFaceIndices(G,'North');
    finx_eb = boundaryFaceIndices(G,'East');
    [n_sb, ~] = gridFaceNodes(G, finx_sb);
    [n_nb, ~] = gridFaceNodes(G, finx_nb);
    [n_eb, ~] = gridFaceNodes(G, finx_eb);
    G.nodes.coords(n_sb,3) = G.nodes.coords(n_sb,3) - opt.dz_s; %0;
    G.nodes.coords(n_nb,3) = G.nodes.coords(n_nb,3) - opt.dz_n; %0;
    G.nodes.coords(n_eb,3) = G.nodes.coords(n_eb,3) - opt.dz_e; %1.05;
    G = computeGeometry(G); % need to recompute geometry
    Gt = topSurfaceGrid(G);

    
    % uncomment the following to plot the grid
%     figure; set(gcf,'Position',[2817 911 560 420])
%     Gt_zshifted = Gt; 
%     Gt_zshifted.nodes.z = Gt_zshifted.nodes.z - 15;
%     plot_opts = {'edgeColor', 'k', 'edgeAlpha', 0.1};
%     plotGrid(G, plot_opts{:});
%     plotCellData(Gt_zshifted, Gt_zshifted.cells.z, plot_opts{:});
%     view(30,25); axis tight
    
    
    % uncomment the following to inspect trapping structure
%     figure; set(gcf,'Position',[3403 846 1517 508])
%     % Cell-based: to be used when comparing trap heights to co2 heights.
%     subplot(1,2,1); plot(1);
%     ta_cell = trapAnalysis(Gt, true);
%     mapPlot(gcf, Gt, 'traps', ta_cell.traps, ...
%         'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
%         'rivers', ta_cell.cell_lines, 'rivercolor', [1 0 0], 'maplines', 20);
%     colorizeCatchmentRegions(Gt, ta_cell);
%     axis equal
%     title('cell-based')
%     % Node-based:
%     subplot(1,2,2); plot(1);
%     ta_node = trapAnalysis(Gt, false);
%     mapPlot(gcf, Gt, 'traps', ta_node.traps, ...
%         'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
%         'rivers', ta_node.cell_lines, 'rivercolor', [1 0 0], 'maplines', 20);
%     colorizeCatchmentRegions(Gt, ta_node);
%     axis equal
%     title('node-based')
%     
%     figure; set(gcf,'Position',[2895 170 1995 605])
%     subplot(1,2,1)
%     plotCellData(Gt, ta_cell.trap_regions); title('cell-based')
%     subplot(1,2,2)
%     plotCellData(Gt, ta_node.trap_regions); title('node-based')


end

