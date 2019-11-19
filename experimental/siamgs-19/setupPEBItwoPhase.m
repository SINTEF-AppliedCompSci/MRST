function [G, rock, fluid, state0, schedule] = setupPEBItwoPhase(varargin)
    
    opt = struct('n'      , 10    , ...
                 'dt'     , 30*day, ...
                 'time'   , 6*year, ...
                 'nRampup', 0     );
    opt = merge_options(opt, varargin{:});
             
    l = 400*meter;
    rng(1);
    d = 1/opt.n;
%     [x,y] = ndgrid(linspace(d/2,1-d/2,opt.n));
%     x     = ([x(:), y(:)] + (2*rand(opt.n^2,2) - 1)*0.25*d )*l;
%     bnd   = [0, 0; 1, 0; 1, 1; 0, 1; 0,0]*l;
%     G     = clippedPebi2D(x, bnd);
%     G     = removeShortEdges(G, 1);
    if 0
    [x,y] = ndgrid(linspace(d/2,1-d/2,opt.n));
    x     = ([x(:), y(:)] + (2*rand(opt.n^2,2) - 1)*0.25*d )*l;
    bnd   = [0, 0; 1, 0; 1, 1; 0, 1; 0,0]*l;
    G     = clippedPebi2D(x, bnd);
    G     = removeShortEdges(G, 1);
    else
    G = pebiGrid(l/opt.n, [l,l]);
    end
    
    
    nLayers = 5;
    layerThickness = diff(linspace(0,l*0.2,nLayers+1));
    G = makeLayeredGrid(G, layerThickness);
    
    G = cartGrid([10,10,3], [l,l,l]);
    
    G = computeGeometry(G);
    G = computeCellDimensions2(G);
    [G.cells.equal, G.faces.equal] = deal(false);
    rock = makeRock(G, 100*milli*darcy, 0.4);
    
    fluid = initSimpleADIFluid('phases', 'WO'                            , ...
                               'rho'   , [1000, 500]*kilogram/meter^3, ...
                               'mu'    , [1   , 0.5]*centi*poise       , ...
                               'n'     , [2   , 2]                   , ...
                               'c'     , 0*[1e-6, 1e-5]/barsa         );
    
    rate   = 0.1*sum(poreVolume(G, rock))/year;

    W = [];
    xc = G.cells.centroids(:,1:2);
    xw = [0,0; l,l];
    dst = pdist2(xc, xw);
    d = mean(G.cells.volumes)^(1/3);
    c = dst(:,1) < d;
    W = addWell(W, G, rock, c, ...
                     'type'  , 'rate', ...
                     'val'   , rate  , ...
                     'comp_i', [0,1] );
    dst = pdist2(xc, xw);
    c = dst(:,2) < d;
    W = addWell(W, G, rock, c, ...
                'type'  , 'bhp'   , ...
                'val'   , 50*barsa, ...
                'comp_i', [0,1]   );
    
    dt = rampupTimesteps2(opt.time, opt.dt);
    schedule = simpleSchedule(dt, 'W', W);
    
    % Initial state
    sW          = 1.0;
    state0      = initResSol(G, 100*barsa, [sW,1-sW]);
    [state0.rs, state0.rv] = deal(0);

end