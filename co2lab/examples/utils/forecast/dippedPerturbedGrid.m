function Gt = dippedPerturbedGrid(varargin)
% Construct grid similar to one used in trappingExample1.m, with optional
% modifications such that trap edges to not coincide with boundary edges.
%
% nx, ny - grid discretization
% Lx, Ly, H - dimensions of grid (in meters)
% dz_s, dz_n, dz_e - amount to chop from south, north, east, etc side

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    [opt.Lx, opt.Ly, opt.H]         = deal(10000, 5000, 50);
    [opt.nx, opt.ny]                = deal(200, 100);
    [opt.dz_s, opt.dz_n, opt.dz_e]  = deal(0, 0, 0);
    opt.originalGrid                = true; % use original grid
    opt = merge_options(opt, varargin{:});
    
    G = cartGrid([opt.nx opt.ny 1],[opt.Lx opt.Ly opt.H]);
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
        % BETA impacts amount of squeeze/stretch in waves
        % PHI impacts curvature of top from north to south
        ALPHA = 0.2;
        WEDGE = 0.15;
        BEND = 0.075;
        SLOPE = 1;
        gamma1 = 6;
        gamma2 = 6;
        beta = 1;
        phi = 1.25;

        fac = 0;
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
    
    
    % Perform some removing of cells
    % (or extract subgrid)
    x = G.cells.centroids(:,1);
    y = G.cells.centroids(:,2);
    xcinxE = find(x >= 0.85*opt.Lx); xcinxW = find(x <= 0.14*opt.Lx);
    ycinxN = find(y >= 0.9*opt.Ly); ycinxS = find(y <= 0.09*opt.Ly);
    cells2cut = unique([xcinxE; xcinxW; ycinxN; ycinxS]);
    %figure; plotGrid(G, 'facecolor','none'); plotCellData(G, ones(G.parent.cells.num,1), cells2cut);
    %G = removeCells(G, cells2cut);
    tmp = ones(G.cells.num,1);
	tmp(cells2cut) = 0;
	cells2keep = find(tmp == 1);
    G = extractSubgrid(G, cells2keep);
    
    
    % Shift G.nodes.coords(x,y) so that grid starts at (0,0)
    G.nodes.coords(:,1) = G.nodes.coords(:,1) - min(G.nodes.coords(:,1));
    G.nodes.coords(:,2) = G.nodes.coords(:,2) - min(G.nodes.coords(:,2));
    
    
    % Alter the external boundaries to ensure trap regions do not touch
    % edges (i.e., forces co2 to spill out along these sides)
    finx_sb = searchForBoundaryFaces(G,'South');
    finx_nb = searchForBoundaryFaces(G,'North');
    finx_eb = searchForBoundaryFaces(G,'East');
    [n_sb, ~] = gridFaceNodes(G, finx_sb);
    [n_nb, ~] = gridFaceNodes(G, finx_nb);
    [n_eb, ~] = gridFaceNodes(G, finx_eb);
    G.nodes.coords(n_sb,3) = G.nodes.coords(n_sb,3) - opt.dz_s; %0;
    G.nodes.coords(n_nb,3) = G.nodes.coords(n_nb,3) - opt.dz_n; %0;
    G.nodes.coords(n_eb,3) = G.nodes.coords(n_eb,3) - opt.dz_e; %1.05;
    
    
    % Recompute geometry, get top surface grid
    G = computeGeometry(G);
    Gt = topSurfaceGrid(G);
    %figure; plotGrid(Gt); view(3);


    
    % uncomment the following to plot the grid
    %{
    figure; set(gcf,'Position',[2817 911 560 420])
    Gt_zshifted = Gt; 
    Gt_zshifted.nodes.z = Gt_zshifted.nodes.z - 15;
    plot_opts = {'edgeColor', 'k', 'edgeAlpha', 0.1};
    plotGrid(G, plot_opts{:});
    plotCellData(Gt_zshifted, Gt_zshifted.cells.z, plot_opts{:});
    view(30,25); axis tight
    %}
    
    %{
    % uncomment the following to inspect trapping structure
    figure; set(gcf,'Position',[3403 846 1517 508])
    % Cell-based: to be used when comparing trap heights to co2 heights.
    subplot(1,2,1); plot(1);
    ta_cell = trapAnalysis(Gt, true);
    mapPlot(gcf, Gt, 'traps', ta_cell.traps, ...
        'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
        'rivers', ta_cell.cell_lines, 'rivercolor', [1 0 0], 'maplines', 20);
    colorizeCatchmentRegions(Gt, ta_cell);
    axis equal
    title('cell-based')
    % Node-based:
    subplot(1,2,2); plot(1);
    ta_node = trapAnalysis(Gt, false);
    mapPlot(gcf, Gt, 'traps', ta_node.traps, ...
        'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
        'rivers', ta_node.cell_lines, 'rivercolor', [1 0 0], 'maplines', 20);
    colorizeCatchmentRegions(Gt, ta_node);
    axis equal
    title('node-based')
    
    figure; set(gcf,'Position',[2895 170 1995 605])
    subplot(1,2,1)
    plotCellData(Gt, ta_cell.trap_regions); title('cell-based')
    subplot(1,2,2)
    plotCellData(Gt, ta_node.trap_regions); title('node-based')
    %}

end
