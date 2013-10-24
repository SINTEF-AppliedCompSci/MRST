mrstModule add coarsegrid gridtools mrst-gui

%% Extract a subset of the fine scale Utsira formation 
grdecl = getAtlasGrid('utsirafm', 'coarsening', 1);
G = processGRDECL(grdecl{1});
G = computeGeometry(G(1));
cdims = G.cartDims;
G = extractSubgrid(G, find(G.cells.centroids(:, 2) > 6.62e6));
G.cartDims = cdims;
G = computeGeometry(G);

Gt = topSurfaceGrid(G);
rock = struct('poro', repmat(.3, G.cells.num, 1), ...
              'perm', repmat(100*milli*darcy, G.cells.num, 1));
%% Compute trapping analysis and find the sorted injection trees
res = trapAnalysis(Gt, 'cell');
trees = maximizeTrapping(Gt, 'res', res);

%% Find the ideal injection spot for each trap region
% Use the largest z value to find the ideal injection spot: This
% corresponds to the deepest possible place in the reservoir surface to
% inject.
N = 5;
W = [];
injRate = sum(G.cells.volumes)/(10000*N*year);

for i = 1:N
    region = find(res.trap_regions == trees(i).root);
    [~, minpoint] = max(Gt.cells.z(region));
    W = addWell(W, G, rock, region(minpoint), ...
        'Name', num2str(i), ...
        'val', injRate, ...
        'Type', 'rate', ...
        'comp_i', [1 0]);
    
end

clf;
v = trapPathsValue(Gt, trees, res);
plotCellData(Gt, log10(v), v>0)
plotGrid(Gt, 'facec', 'none', 'edgea', .05)
plotWell(G, W)
axis tight off
view(115, 60)




%% Set up fluid
muw = 0.30860;  rhow = 975.86; sw    = 0.1;
muc = 0.056641; rhoc = 686.54; srco2 = 0.2;
kwm = [0.2142 0.85];

mu  = [muc  muw ] .* centi*poise;
rho = [rhoc rhow] .* kilogram/meter^3;

fluid = initSimpleVEFluid_s('mu' , mu , 'rho', rho, ...
                                'height'  , Gt.cells.H,...
                                'sr', [srco2, sw],'kwm',kwm);
fluid.sr = srco2;
fluid.sw = sw;

%% Define vertical pressure distribution
topPressure = 350*barsa;
gravity on;
grav     = gravity();
topPos   = min(Gt.cells.z);
pressure = @(z) topPressure + rho(2)*(z - topPos)*grav(3);

%% Schedule
T_injection = 100*year;
T_migration = 500*year;
Ni = 20;
Nm = 100;

T_tot = T_injection + T_migration;
dTi   = T_injection / Ni;  % short time steps during injection
dTm   = T_migration / Nm;  % longer steps during migration

rock2D    = averageRock(rock, Gt);
T = computeTrans(Gt, rock2D);
T = T.*Gt.cells.H(gridCellNo(Gt));


% Add pressure boundary 
bnd = boundaryFaces(Gt);
bc = addBC([], bnd, 'pressure', pressure(Gt.faces.z(bnd)), 'sat', [0 1]);

% Convert to 2D wells
W2D = convertwellsVE_s(W, G, Gt, rock2D,'ip_tpf');

%%  Set up initial reservoir conditions
% The initial pressure is set to hydrostatic pressure. Setup and plot.
sol = initResSolVE_s(Gt, pressure(Gt.cells.z), 0);
sol.wellSol = initWellSol(W2D, 0);
sol.h = zeros(Gt.cells.num, 1);


opts = {'Saxis',     [0 1-fluid.sw], ...
        'plotPlume', true, ...
        'wireH', true, ...
        'wireS', true, ...
        'slice',     [150 800], ...
        'maxH',      (max(Gt.cells.z) - min(Gt.cells.z))/3, ...
        'plotHist',  true};


fh = plotPanelVE(G, Gt, W, sol, 0, ...
   [volumesVE(Gt, sol, rock2D, fluid, res) 0], opts{:});


%% Run the simulation
W2D = convertwellsVE_s(W, G, Gt, rock2D,'ip_tpf');
% Solve for all timesteps, and plot the results at each timestep.
t = 0;
tt = ' (Injecting)';
totVol = 0;

% Estimate total timesteps
tstep_tot = ceil(T_injection/dTi) + ceil(T_migration/dTm); 
i = 1;

% Waitbar to show progress
hwbar = waitbar(0, 'Initializing...');
wbar = @(i, t, status) waitbar(t/T_tot, hwbar, ...
   sprintf('Timestep %d of %d, T=%s%s', i, tstep_tot, ...
   formatTimeRange(floor(t/year)*year), status));
while t < T_tot;
    if ishandle(hwbar)
        wbar(i, t, tt);
    end
    if t >= T_injection 
        W2D = [];
        dT  = dTm;
        % Do a shorter timestep
        dT = min(dT, T_tot - t);
        tt  = ' (Migrating)';
    else
        dT = dTi;
        dT = min(dT, T_injection - t);
    end
    sol = incompTPFA(sol, Gt, T, fluid, 'wells', W2D, 'bc', bc);
    sol = implicitTransport(sol, Gt, dT, rock2D, fluid, ...
                            'wells', W2D, 'bc', bc, 'Verbose', false);
    t = t + dT;

    % Plotting
    [s h hm] = normalizeValuesVE(Gt, sol, fluid);
    sol.h = h;
    sol.h_max = hm;

    if true
        % Use advanced plotting
        if ~isempty(W2D)
            totVol = totVol + sum([W2D.val])*dT;
        end
        vol = volumesVE(Gt, sol, rock2D, fluid, res);
        plotPanelVE(G, Gt, W, sol, t, [vol totVol], opts{:}); 
    else
        set(0,'CurrentFigure', fh);
        [a,b] = view();
        clf
        plotCellData(G, s, 'edgec', 'k', 'edgea', .1, 'edgec', [.6 .6 .6]);
        plotWell(G, W); caxis([0 .9]);
        title([formatTimeRange(t) tt])
        colorbar
        axis tight off
        view(a,b);
        drawnow
    end
    i = i + 1;
end
