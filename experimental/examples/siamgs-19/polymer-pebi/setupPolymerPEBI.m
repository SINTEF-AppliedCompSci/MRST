function [G, rock, state0, schedule] = setupPolymerPEBI(fluid, varargin)
    
    opt = struct('n'      , 10    , ...
                 'dt'     , 30*day, ...
                 'nRampup', 0     );
    opt = merge_options(opt, varargin{:});
             
    l = 400*meter;
    rng(1);
    d = 1/opt.n;
    [x,y] = ndgrid(linspace(d/2,1-d/2,opt.n));
    x     = ([x(:), y(:)] + (2*rand(opt.n^2,2) - 1)*0.25*d )*l;
    bnd   = [0, 0; 1, 0; 1, 1; 0, 1; 0,0]*l;
    G     = clippedPebi2D(x, bnd);
%     nLayers = 5;
%     layerThickness = diff(linspace(0,l*0.2,nLayers+1));
%     G = makeLayeredGrid(G, layerThickness);
    G = computeGeometry(G);
    G = computeCellDimensions2(G);
    [G.cells.equal, G.faces.equal] = deal(false);
    rock = makeRock(G, 100*milli*darcy, 0.4);
             
    opt = merge_options(opt, varargin{:});
    
    wTime1 = 0.5*year;
    wTime2 = 10*year;
    pTime  = 0.5*year;
    rate   = 0.1*sum(poreVolume(G, rock))/year;

    W = [];
    xc = G.cells.centroids(:,1:2);
    xw = [0,0; l,l];
    dst = pdist2(xc, xw);
    d = mean(G.cells.volumes)^(1/2);
    c = dst(:,1) < d;
    W = addWell(W, G, rock, c, ...
                     'type'  , 'rate', ...
                     'val'   , rate  , ...
                     'comp_i', [1,0] );
    dst = pdist2(xc, xw);
    c = dst(:,2) < d;
    W = addWell(W, G, rock, c, ...
                'type'  , 'bhp'   , ...
                'val'   , 50*barsa, ...
                'comp_i', [1,0]   );
    [W.c]    = deal(0);
    schedule = makePolymerSlugSchedule(W, fluid, ...
                                'dt'     , opt.dt          , ...
                                'pTime'  , pTime           , ...
                                'wTime'  , [wTime1, wTime2], ...
                                'nRampup', opt.nRampup     );
    % Initial state
    sW          = 0.18;
    state0      = initResSol(G, 100*barsa, [sW,1-sW]);
    [state0.c, state0.cmax] = deal(zeros([G.cells.num, 1]));

end