clear
close all
mrstModule add ad-core ad-blackoil ad-props mrst-gui wellpaths deckformat
%%
GC = cartGrid([25, 25], [200, 200]);
GC = computeGeometry(GC);

ns = 12;
th = linspace(0.85, 0.7, ns+1)' * pi;
traj = [150*cos(th)+200, 150*sin(th)];

pbdy = [   136   150
    145    95
    90    30
    50    50
    45   105
   90   160];

cCtro = GC.cells.centroids;
in = inpolygon(cCtro(:,1), cCtro(:,2), pbdy(:,1), pbdy(:,2));
CI = find( in );
CO = find(~in );

[bnv, bcv] = demoGetBdyNodesCells(GC, CI);

figure, hold on, axis equal tight off
plotGrid(GC, CO, 'facecolor', 'none')
plotGrid(GC, CI, 'facecolor', 'y')
demoPlotLine(traj, 'ko-', 'b', 4)
demoPlotPoly(pbdy, 'k^-', 'r', 5)
demoPlotPoly(GC.nodes.coords(bnv,:), 'ks-', 'g', 4)
legend('GC outside the VOI', 'GC inside the VOI', 'Well path' ,...
    'Specified VOI boundary','Clipped VOI boundary')
%%
ly = 12;
ny = 8;
na = 6;
pw0 = arrayfun(@(ii)pointsSingleWellNode(traj, ly, ny, na, ii), (1:ns+1)');
pw  = [vertcat(pw0.cart); vertcat(pw0.rad)];
assert(all(inpolygon(pw(:,1), pw(:,2), pbdy(:,1), pbdy(:,2))),...
    ['Points outside the boundary were detected, please reduce ',...
    'the size of Cartesian region']);
[tw, ~, bnw] = getConnListAndBdyNodeWR2D(pw0, ny, na);

GWR = tessellationGrid(pw, tw);

figure, hold on, axis equal tight off
plotGrid(GC, CO, 'facecolor', 'none')
plotGrid(GWR, 'facecolor', 'c')
demoPlotPoly(GC.nodes.coords(bnv, :), 'ks-', 'g', 4)
demoPlotPoly(pw(bnw, :), 'ko-', 'y', 3)
legend('GC outside the VOI', 'GWR (WR grid)', 'VOI boundary' ,'WR boundary')
%%
pob  = GC.nodes.coords(bnv,:);
pob2 = GC.cells.centroids(bcv, :);
pib  = pw(bnw, :);
[pIn, pOut, R] = demoComputeAuxPts(pw, bnw, 0.35);
pib2 = pOut;

figure, hold on, axis equal tight off
demoPlotPoly(pob,  'ks-', 'g', 4)
demoPlotPoly(pib,  'ko-', 'y', 3)
demoPlotPoly(pob2, 'ks-', 'r', 4)
demoPlotPoly(pib2, 'ko-', 'b', 3)
legend('VOI boundary' ,'WR boundary', 'VOI boundary for tri. pts. generation' ...
    ,'WR boundary for tri. pts. generation')
%%
% Generate basic points from distmesh
[pdis, tdis] = passToDistmesh(pib2, pob2, 0.2, 500, 'pIBRadius', R);

pauxw = pIn;
pauxv = GC.cells.centroids(CO, :);

% Get Voronoi points and connectivity list
pall = [pdis; pauxv; pauxw];
[pVor, tVor] = voronoin(pall, {'Qbb','Qz'});
%%
figure, hold on, axis equal off
voronoi(pall(:,1), pall(:,2), '-')
demoPlotPoly(pdis,  'ro',  'r', 2)
demoPlotPoly(pauxv, 'bo',  'b', 2)
demoPlotPoly(pauxw, 'k^',  'k', 3)
demoPlotPoly(pob,   'ks-', 'g', 4)
demoPlotPoly(pib,   'ko-', 'y', 3)
legend('', 'Initial Voronoi digram', 'Basic sites' ,'Auxiliary sites (CO)', ...
    'Auxiliary sites (WR)',  'VOI boundary' ,'WR boundary')
xlim([-20, 220])
ylim([-20, 220])
%%
[pv, tv, bnv2] = demoHandleVoronoiDiag(pVor, tVor, pib, pob, pw, tw, bnw);

tv = sortPtsCounterClockWise(pv, tv);
GV = tessellationGrid(pv, tv);
GV = computeGeometry(GV);

figure, hold on, axis equal tight off
plotGrid(GC, CO, 'facecolor', [.0, .44, .74])
plotGrid(GV, 'facecolor', [.85, .32, .09])
legend('GC outside the VOI', 'GV (New VOI grid)')
%%
[GC, ~, ~, mapn] = removeCells(GC, CI);
%%
bnv = arrayfun(@(n)find(mapn == n), bnv);
%%
nNo = (1:GC.nodes.num)';
idx = ~ismember(nNo, bnv);
resn = (1:nnz(idx))' + size(pv,1);
nNo(idx)  = resn;
nNo(bnv)  = bnv2;
%%
pc = GC.nodes.coords;
pc = pc(idx, :);
%%
[cnodes, pos] = gridCellNodes(GC, (1:GC.cells.num));
cnodes = nNo(cnodes);
tc = arrayfun(@(c)cnodes(pos(c):pos(c+1)-1), (1:GC.cells.num)', 'UniformOutput', false);
%%
p = [pv; pc];
t = [tv; tc];
%%
t = sortPtsCounterClockWise(p, t);
G = tessellationGrid(p, t);
G = computeGeometry(G);
%%
figure, hold on
plotGrid(G)
plotFaces(G, ~all(G.faces.neighbors, 2), 'edgecolor', 'r')
%%
G2 = assembleGrids({GV, GC});
figure, hold on
plotGrid(G2)
plotFaces(G2, ~all(G2.faces.neighbors, 2), 'edgecolor', 'r')
%%
bfnc = [bnv, bnv([2:end, 1])] + GV.nodes.num;
bfnv = [bnv2, bnv2([2:end, 1])];

[fn, pos] = gridFaceNodes(G2, (1:G2.faces.num)');
assert(all(diff(pos)==2))
fn = reshape(fn, 2, [])';
%%
bfnc = sort(bfnc, 2);
bfnv = sort(bfnv, 2);
fn = sort(fn, 2);
%%
ismember(fn, bfnc, 'rows')
%%
bfc = arrayfun(@(i)find(ismember(fn, bfnc(i,:), 'rows')), (1:size(bfnc,1))');
bfv = arrayfun(@(i)find(ismember(fn, bfnv(i,:), 'rows')), (1:size(bfnv,1))');
assert(all(~all(G2.faces.neighbors([bfc; bfv], :),2)))
%%
f_nnc = [bfc, bfv];
c_nnc = [sum(G2.faces.neighbors(bfc,:),2), sum(G2.faces.neighbors(bfv,:),2)];

figure, hold on, axis equal tight off
plotGrid(G2, 'facecolor', 'none')
arrayfun(@(i)demoPlotLine(G2.cells.centroids(c_nnc(i, :), :)...
    , 'ks-', 'g', 4), (1:size(c_nnc,1))');
plotFaces(G2, ~all(G2.faces.neighbors, 2), 'edgecolor', 'r', 'linewidth', 1)
%%
rock.perm = 100*ones(G.cells.num,2) * (milli*darcy);
cWR = (1:numel(tw))';
rock.perm(cWR, :) = rock.perm(cWR, :) * 10;
rock.poro = 0.2*ones(G.cells.num,1);
%%
figure
subplot(2,2,1), axis equal tight off
plotCellData(G, (1:G.cells.num)')
title('Cell indices of G')
subplot(2,2,2), axis equal tight off
plotCellData(G2, (1:G2.cells.num)')
title('Cell indices of G2')
subplot(2,2,3), axis equal tight off
plotCellData(G, rock.perm(:,1))
title('X-Perm. of G')
subplot(2,2,4), axis equal tight off
plotCellData(G2, rock.perm(:,1))
title('X-Perm. of G2')
%%
fluid = initSimpleADIFluid('mu',    [1, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0], ...
                           'c',     [1e-5, 1e-4, 0] * (barsa)^(-1), ...
                           'cR',    1e-5 * (barsa)^(-1));
                       
gravity reset on
model  = TwoPhaseOilWaterModel(G,  rock, fluid);
model2 = TwoPhaseOilWaterModel(G2, rock, fluid);

N = G2.faces.neighbors;
intInd = all(N,2);
N = N(intInd, :);
N = [N; c_nnc];
T_all = model2.operators.T_all;
T = T_all(intInd);
hT_nnc = T_all(f_nnc);
T_nnc = 1./(1./hT_nnc(:,1) + 1./hT_nnc(:,2));
T = [T; T_nnc];
model2.operators = setupOperatorsTPFA(G2, rock, 'neighbors', N, 'trans', T);
model2.operators.T_all = [T_all; T_nnc];

%% Define initial state
sW = zeros(G.cells.num, 1);
sat = [sW, 1 - sW];
state0 = initResSol(G, 200*barsa, sat);
%% Define wells and simulation schedule
simTime = 2*year;
nstep   = 20;
refine  = 10;

% Producer

W = addWell([], G, rock, G.cells.num, 'Name', 'PROD', 'sign', -1, ...
    'comp_i', [1, 1], 'Val', 150*barsa, 'Type', 'bhp');

% Injector
pv      = poreVolume(G, rock);
injRate = 0.5*sum(pv)/simTime;
D = sqrt(sum(G.cells.centroids.^2, 2));
wc = find( D==min(D) );
W = addWell(W, G, rock, wc, 'Name', 'INJ', 'sign', 1, ...
    'comp_i', [1, 0], 'Val', injRate, 'Type', 'rate');      

% Define the timesteps
% Compute the timesteps
startSteps = repmat((simTime/(nstep + 1))/refine, refine, 1);
restSteps =  repmat(simTime/(nstep + 1), nstep, 1);
timesteps = [startSteps; restSteps];

% Set up the schedule containing both the wells and the timesteps
schedule = simpleSchedule(timesteps, 'W', W);
%%
figure, hold on
plotGrid(G, 'facecolor', 'none')
demoPlotLine(G.cells.centroids(W(1).cells,:), 'ko', 'r', 6)
demoPlotLine(G.cells.centroids(W(2).cells,:), 'k^', 'b', 6)
%%
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule);
[wellSols2, states2, report2] = simulateScheduleAD(state0, model2, schedule);
%%
plotWellSols({wellSols, wellSols2}, report.ReservoirTime)
%%
t = 30;
figure
subplot(2,2,1), axis equal tight off
plotCellData(G, states{t}.pressure/barsa, 'edgecolor', 'none')
subplot(2,2,2), axis equal tight off
plotCellData(G2, states2{t}.pressure/barsa, 'edgecolor', 'none')
title('Cell indices of G')
subplot(2,2,3), axis equal tight off
plotCellData(G, states{t}.s(:,2), 'edgecolor', 'none')
subplot(2,2,4), axis equal tight off
plotCellData(G2, states2{t}.s(:,2), 'edgecolor', 'none')
title('Cell indices of G')
