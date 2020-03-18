
%% run field model
clear all; close all; clc;

mrstModule add incomp mimetic ad-blackoil ad-core ad-props mpfa nfvm ...
    blackoil-sequential wellpaths mrst-gui

nx = 7;
ny = 5;
nz = 3;
refine = 2;
G = cartGrid([nx, ny, nz]*refine, [1, 2, 3]);
%G = tetrahedralGrid(G.nodes.coords);
%G = twister(G, 0.1);
G = computeGeometry(G);


% %% 2D
% dims = [21, 10];
% G = cartGrid(dims, [2, 1]);
% makeSkew = @(c) c(:, 1) + .4 * (1 - (c(:, 1) - 1).^2) .* (1 - c(:, 2));
% G.nodes.coords(:, 1) = 2 * makeSkew(G.nodes.coords);
% G.nodes.coords(:, 1) = G.nodes.coords(:, 1) * 1000;
% G.nodes.coords(:, 2) = G.nodes.coords(:, 2) * 1000;
%
% %G = twister(G, 0.1);
% G = computeGeometry(G);

% G = squareGrid([4, 4], [3, 3], 'grid_type', 'mixed3', 'disturb', 0.01);
% G = computeGeometry(G);


% % Extract faces in x direction for plotting fluxes
% tol = 1e-5;
% pts = [tol, 0.5 + tol, 0.5 + tol; ...
%     1 - tol, 0.5 + tol, 0.5 + tol];
% wellpath = makeSingleWellpath(pts);
% %plotWellPath(wellpath);
% cells_xdir = findWellPathCells(G, wellpath);
% [Gsub, ~, subgf] = extractSubgrid(G, cells_xdir);

rock = makeRock(G, 1, 1);
pv = sum(poreVolume(G, rock));

mu_value = 1;
rho_value = 1;
fluid = initSingleFluid('mu', mu_value, 'rho', rho_value);
fluidad = initSimpleADIFluid('phases', 'W', 'mu', mu_value, 'rho', rho_value);
fluid_glpk = initSimpleADIFluid();
p0 = 1.0;
s0 = 1.0;

% Wells
W = [];
% cellsWell1 = well_cells(G, '1');
% cellsWell2 = well_cells(G, '2');

% radius = 0.05/10;
% W = addWell(W, G, rock, cellsWell1, 'Type', 'rate', ...
%             'InnerProduct', 'ip_tpf', ...
%             'comp_i', [1], ...
%             'Val', 1.0/day(), 'Radius', radius, 'Name', 'I');
% W = addWell(W, G, rock, cellsWell2, 'Type', 'bhp', ...
%             'InnerProduct', 'ip_tpf', ...
%             'comp_i', [1], ...
%             'Val', 1.0e5, 'Radius', radius, 'Dir', 'y', 'name', 'P');

% W = addWell(W, G, rock, cellsWell1, 'Type', 'bhp', ...
%             'InnerProduct', 'ip_tpf', ...
%             'comp_i', [1], ...
%             'Val', -1.0e5, 'Radius', radius, 'Dir', 'y', 'name', 'P');

% % Symmetric well pattern
% [ii, jj] = gridLogicalIndices(G);
% % Injector + two producers
% W = [];
% W = addWell(W, G, rock, find(ii == ceil(G.cartDims(1)/2) & jj == G.cartDims(2)), 'comp_i', [1], 'type', 'rate', 'val', pv/year);
% W = addWell(W, G, rock, find(ii == G.cartDims(1) & jj == 1), 'comp_i', [1], 'type', 'bhp', 'val', 50*barsa);
% W = addWell(W, G, rock, find(ii == 1 & jj == 1), 'comp_i', [1], 'type', 'bhp', 'val', 50*barsa);


bc_flow_ntpfa.face = boundaryFaces(G);
bc_flow_ntpfa.type = repmat({'flux'}, numel(bc_flow_ntpfa.face), 1);
bc_flow_ntpfa.value = repmat({@(x) 0}, numel(bc_flow_ntpfa.face), 1);

tol = 1e-5;
xmin = min(G.nodes.coords(:, 1));
ix = find(G.faces.centroids(bc_flow_ntpfa.face, 1) < xmin+tol);
bc_flow_ntpfa.type(ix) = repmat({'pressure'}, numel(ix), 1);
bc_flow_ntpfa.value(ix) = repmat({@(x) 1}, numel(ix), 1);

xmax = max(G.nodes.coords(:, 1));
ix = find(G.faces.centroids(bc_flow_ntpfa.face, 1) > xmax-tol);
bc_flow_ntpfa.type(ix) = repmat({'pressure'}, numel(ix), 1);
bc_flow_ntpfa.value(ix) = repmat({@(x) 2}, numel(ix), 1);

bc_std = [];
if isfield(G, 'cartDims')
    bc_std = pside(bc_std, G, 'xmin', 1);
    bc_std = pside(bc_std, G, 'xmax', 2);
else
    disp('No bc set')
end

% bc_flow_ntpfa.face = boundaryFaces(G);
% bc_flow_ntpfa.type = repmat({'flux'}, [numel(bc_flow_ntpfa.face), 1]);
% bc_flow_ntpfa.value = repmat({@(x)0}, [numel(bc_flow_ntpfa.face), 1]);
% bc_std = [];

state0 = initState(G, W, p0, s0);
state0_glpk = initResSol(G, p0, [0, s0]);

dt = 1;
schedule = simpleSchedule(dt, 'W', W, 'bc', bc_std);

results = {};

model = GenericBlackOilModel(G, rock, fluidad, ...
    'water', true, 'oil', false, 'gas', false);

%% MPFA Xavier
% solvers{end+1} = 'mpfa';
% disp(solvers{end});
% %mpfastruct = computeNeumannMultiPointTrans(G, rock);
% %state = incompMPFA3(G, mpfastruct, W, 'outputFlux', true);
% mpfastruct = computeMultiPointTrans2(G, rock);
% state = incompMPFAbc(G, mpfastruct, bc_std, 'outputFlux', true);
% figure,plotCellData(G,state.pressure/1e6);plotWell(G,W);view(G.griddim)
% title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';
% results{end+1} = state;

%% TPFA from NFVM
% solvers{end+1} = 'tpfa';
% disp(solvers{end});
% state=FlowTPFA(G,TransTPFA(G,rock,bc_flow_ntpfa),fluid,W);
% figure,plotCellData(G,state.pressure/1e6);plotWell(G,W);view(G.griddim)
% title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';
% Qinj=sum(state.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
% state=tracerTransport_implicit(G,rock,W,state,dt,nstep);
% figure, plotCellData(G,state.cres(:,end));plotWell(G,W);
% view(G.griddim);colorbar;title(['Concentration ', solvers{end}]);
% drawnow
% results{end+1} = state;

%% incomp TPFA
% solvers{end+1} = 'incomp tpfa';
% disp(solvers{end});
% T = computeTrans(G,rock);
% state = incompTPFA(state0, G, T, fluid, 'wells', W, 'bc', bc_std, ...
%                    'matrixOutput', true);
% figure,plotToolbar(G,state);plotWell(G,W);view(G.griddim)
% title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';
% results{end+1} = state;

%% TPFA new framework
results{end+1}.tag = 'MRST TPFA new framework';
fprintf('\n%s\n', results{end}.tag);
[wellSols, states] = simulateScheduleAD(state0, model, schedule);
figure, plotToolbar(G, states);
%plotWell(G, W);
view(G.griddim)
title(['Pressure ', results{end}.tag]);
h = colorbar;
h.Label.String = 'Pressure [Pa]';
results{end}.state = states{1};

%% NTPFA Kobaisi in new framework
results{end+1}.tag = 'NTPFA new framework';
fprintf('\n%s\n', results{end}.tag);
mvstatus = mrstVerbose;
mrstVerbose on
model_ntpfa = setNTPFADiscretization(model);
[wellSols, states] = simulateScheduleAD(state0, model_ntpfa, schedule);
figure, plotToolbar(G, states);
%plotWell(G, W);
view(G.griddim)
title(['Pressure ', results{end}.tag]);
h = colorbar;
h.Label.String = 'Pressure [Pa]';
mrstVerbose(mvstatus);
results{end}.state = states{1};

%% NTPFA Kobaisi original
results{end+1}.tag = 'NTPFA Kobaisi';
mvstatus = mrstVerbose;
mrstVerbose on
disp(results{end}.tag);
interpFace = findHAP(G, rock, bc_flow_ntpfa);
interpFace = correctHAP(G, interpFace);
OSflux = findOSflux(G, rock, bc_flow_ntpfa, interpFace);
state = FlowNTPFA(G, bc_flow_ntpfa, fluid, W, OSflux, p0*ones(G.cells.num, 1), 1e-7, 1000, 'matrixOutput', false);
figure, plotToolbar(G, state, 'field', 'pressure');
%plotWell(G, W);
view(G.griddim)
title(['Pressure ', results{end}.tag]);
h = colorbar;
h.Label.String = 'Pressure [Pa]';
results{end}.state = state;
mrstVerbose(mvstatus);

%% NTPFA glpk
% solvers{end+1} = 'NTPFA-glpk';
% disp(solvers{end});
% mrstVerbose on
% mrstDebug on
% model = PressureOilWaterModelNTPFAopt(G,rock,fluid_glpk);
% state = incompSinglePhaseNTPFA(model,state0_glpk,'bc', bc_glpk,'src', []);
% figure,plotToolbar(G,state);plotWell(G,W);view(G.griddim)
% title(['Pressure ', solvers{end}]);h=colorbar;h.Label.String='Pressure [Pa]';
% results{end+1} = state;

%% nonlinear MPFA from Kobaisi
% solvers{end+1} = 'nonlinear MPFA';
% disp(solvers{end});
% %state=NewtonMPFA(G,fluid,W,OSflux,1e-9,100);
% state=FlowNMPFA(G,bc_flow_ntpfa,fluid,W,OSflux,1e-9,10000);
% figure,plotCellData(G,state.pressure/1e6);plotWell(G,W);axis image;view(G.griddim)
% title('Pressure-NMPFA');h=colorbar;h.Label.String='Pressure [Pa]';
% Qinj=sum(state.wellSol(1).flux);tt=1.5*pv/Qinj;nstep=100;dt=tt/nstep;
% state=tracerTransport_implicit(G,rock,W,state,dt,nstep);
% figure, plotCellData(G,state.cres(:,end));plotWell(G,W);
% view(G.griddim);colorbar;title(['Concentration ', solvers{end}]);
% results{end+1} = state;

%% MPFA new framework
results{end+1}.tag = 'MPFA new framework';
model_mpfa = setMPFADiscretization(model);
[wellSols, states] = simulateScheduleAD(state0, model_mpfa, schedule);
figure, plotToolbar(G, states);
%plotWell(G, W);
view(G.griddim)
title(['Pressure ', results{end}.tag]);
h = colorbar;
h.Label.String = 'Pressure [Pa]';
mrstVerbose(mvstatus);
results{end}.state = states{1};

%% MFD (Must be at the end)
results{end+1}.tag = 'MFD';
fprintf('\n%s\n', results{end}.tag);
state = incompMimetic(state0, G, computeMimeticIP(G, rock), fluid, 'wells', W, 'bc', bc_std, 'matrixOutput', false);
%figure, plotCellData(G,state.pressure/1e6);view(G.griddim);
figure, plotToolbar(G, state);
%plotWell(G, W);
view(G.griddim)
title(['Pressure ', results{end}.tag]);
h = colorbar;
h.Label.String = 'Pressure [Pa]';
%figure,plotFaceData(G,cells_xdir, state.flux);colorbar,view(G.griddim)
%title(['flux ', solvers{end}])
results{end}.state = state;

%% primitive comparisons with MFD (must be last)
fprintf('\nCompute some errors\n')
L2err = @(uh, u, vol) num2str(sqrt(sum((uh - u).^2.*vol))); % ./ sqrt(sum(u.^2.*vol)));
for i = 1:numel(results) - 1 % last is mfd
    fprintf('\n%s\n', results{i}.tag);
    figure
    subplot(2, 1, 1)
    diff_pressure = abs(results{i}.state.pressure-results{end}.state.pressure) ./ abs(results{end}.state.pressure);
    plotCellData(G, diff_pressure);
    colorbar;
    view(G.griddim)
    titlename = ['relative error pressure ', results{i}.tag, ' vs ', results{end}.tag];
    title(titlename)
    disp(titlename)
    disp(['diff p max mean L2: ', ...
        num2str(max(diff_pressure)), ' ', ...
        num2str(mean(diff_pressure)), ' ', ...
        L2err(results{i}.state.pressure, results{end}.state.pressure, G.cells.volumes)])

    subplot(2, 1, 2)
    diff_flux = abs(results{i}.state.flux-results{end}.state.flux); %./abs(results{2}.flux);
    plotFaceData(G, diff_flux);
    colorbar;
    view(G.griddim)
    titlename = ['absolute error flux ', results{i}.tag, ' vs ', results{end}.tag];
    title(titlename)
    disp(titlename)
    disp(['diff flux max mean L2: ', ...
        num2str(max(diff_flux)), ' ', ...
        num2str(mean(diff_flux)), ' ', ...
        L2err(results{i}.state.flux, results{end}.state.flux, G.faces.areas)])
end
