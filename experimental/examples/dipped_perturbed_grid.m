function Gt = dipped_perturbed_grid(varargin)
% Construct grid used in /co2lab/examples/papers/CAGEO-75/trappingExample1.m
%
% nx, ny - grid discretization (100x100 or 50x25 etc)
% Lx, Ly, H - dimensions of grid (in meters)
% Grid is 10 km by 5 km by 50 meters deep. 
%

    [opt.Lx, opt.Ly, opt.H] = deal(10000, 5000, 50);
    [opt.nx, opt.ny] = deal(50, 25);
    
    opt = merge_options(opt, varargin{:});
    
    G  = cartGrid([opt.nx opt.ny 1],[opt.Lx opt.Ly opt.H]);
    x  = G.nodes.coords(1:G.nodes.num/2,1)/opt.Lx; 
    y  = G.nodes.coords(1:G.nodes.num/2,2)/opt.Ly;
    z  = G.nodes.coords(1:G.nodes.num/2,3)/opt.H;
    zt = z + x - 0.2*sin(5*pi*x).*sin(5*pi*y.^1.5) - 0.075*sin(1.25*pi*y) + 0.15*sin(x+y);
    % ALPHA impacts amplitude of perturbed waves on top surface
    % WEDGE impacts aspect ratio of domain as a wedge (from north to south)
    % BEND impacts curve of top surface (from east to west)
    % SLOPE impacts top surface's dipping angle
    ALPHA = 0.3;
    WEDGE = 0.5; 
    BEND = 0.1;
    SLOPE = 1.2;
    %zt = z + SLOPE*x - ALPHA*sin(5*pi*x).*sin(5*pi*y.^1.5) - BEND*sin(1.25*pi*y) + WEDGE*sin(x+y);
    zb = 1 + x;
    G.nodes.coords(:,3) = [zt; zb]*opt.H+1000;
    G = computeGeometry(G);
    Gt = topSurfaceGrid(G);

    % uncomment the following to plot the grid
    figure;
    Gt_zshifted = Gt; 
    Gt_zshifted.nodes.z = Gt_zshifted.nodes.z - 15;
    plot_opts = {'edgeColor', 'k', 'edgeAlpha', 0.1};
    plotGrid(G, plot_opts{:});
    plotCellData(Gt_zshifted, Gt_zshifted.cells.z, plot_opts{:});
    view(30,25); axis tight


end

