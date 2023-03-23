function [G, rock, fluid, state0, schedule] = setupThreephasePEBI(varargin)
%Undocumented Utility Function

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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
    
    opt = struct('n'      , 10    , ...
                 'dt'     , 30*day, ...
                 'time'   , 6*year, ...
                 'nRampup', 0     );
    opt = merge_options(opt, varargin{:});
             
    l = 400*meter;
    rng(1);
    d = 1/opt.n;
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
    G = computeGeometry(G);
    G = computeCellDimensions2(G);
    [G.cells.equal, G.faces.equal] = deal(false);
    rock = makeRock(G, 100*milli*darcy, 0.4);
    
    fluid = initSimpleADIFluid('phases', 'WOG'                            , ...
                               'rho'   , [1000, 800, 300]*kilogram/meter^3, ...
                               'mu'    , [0.5 , 1, 0.2]*centi*poise       , ...
                               'n'     , [2   , 2,   2]                   , ...
                               'c'     , [1e-6, 1e-5, 1e-4]/barsa         );
    
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
                     'comp_i', [0,0,1] );
    dst = pdist2(xc, xw);
    c = dst(:,2) < d;
    W = addWell(W, G, rock, c, ...
                'type'  , 'bhp'   , ...
                'val'   , 50*barsa, ...
                'comp_i', [0,0,1]   );
    
    dt = rampupTimesteps2(opt.time, opt.dt);
    schedule = simpleSchedule(dt, 'W', W);
    
    % Initial state
    sW          = 0.18;
    state0      = initResSol(G, 100*barsa, [sW,1-sW, 0]);
    [state0.rs, state0.rv] = deal(0);

end
