%% Example of simulation on a radial-Cartesian hybrid grid
%
% This example demonstrates how to perform simulation on a radial-Cartesian 
% hybrid grid. The additional treatments for such grid involve the 
% calculation of radial transmissibility and the definition of the well
% structure. Also, the simulation performances between the standalone 
% Cartesian grid that uses the Peaceman well model and the hybrid grid are 
% compared.

clear
mrstModule add nwm ad-core ad-blackoil ad-props

%% Build the hybrid grid
% Build the Cartesian grid
GC = cartGrid([15, 15], [1000, 1000]);
GC = computeGeometry(GC);

% Define the well region by logical indices
ij = gridLogicalIndices(GC);
idxI = ij{1} >= 6 & ij{1} <= 10 & ij{2} >= 6 & ij{2} <= 10;
cI = find(idxI); 

% Place the well at the region center
pW  = [500, 500];

% Define radial parameters
[nR, rW, rM] = deal(10, 0.2, 100);

% Get the hybrid grid
G = radCartHybridGrid(GC, cI, rW, rM, nR, pW);
figure, hold on; axis equal off, plotGrid(G)

%% Compute the radial half transmissibility factor
% The radial transmissibility is derived from the radial/angular two-point 
% flux approximation. The derivation assumes the steady-state flow and
% computes 'transmissibility center' by integral of the pressure within the
% area of cell. The half transmissibility = permeability .* factor.

% Assign the radial grid
GR = G.subGrids{1};
% The irregular outermost radial cells are not involved in the calculations
GR.radDims = [GR.radDims(1), GR.radDims(2)-1, 1];

% Compute the factor
skin = 0;
ft = computeRadTransFactor(GR, pW, skin);
ft = ft(~isnan(ft));

%% Assign the radial transmissibility
% Define a homogeneous rock 
rock = makeRock(G, [.25, .25]*darcy, 0.25);

% Linear half transmissibility
hT = computeTrans(G, rock);

% Indices of radial cells involved in the factor calculations
radCells = (1:GR.cells.num-GR.radDims(1))';

% Note the factor is suitable for only isotropic permeability
permR = rock.perm(radCells,1);
n = diff(G.cells.facePos);
permR = rldecode(permR, n(radCells));

% Compute the radial half transmissibility
hTR = permR .* ft;

% Assign the radial half transmissibility
idx = (1:numel(hTR))';
hT(idx) = hTR;

% Compute the full transmissibility
cf = G.cells.faces(:,1);
nf = G.faces.num;
assert( numel(hT) == numel(cf) )
T_all  = 1 ./ accumarray(cf, 1./hT, [nf, 1]);
intCon = all(G.faces.neighbors, 2);
T = T_all(intCon);

%% Setup simulation model
% Define a two phase (oil-water) fluid
fluid = initSimpleADIFluid('mu',    [1, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0],...
                           'c',     [1e-5, 1e-3, 0] * (barsa)^(-1));
                       
% Setup the simulation model                      
model  = TwoPhaseOilWaterModel(G,  rock, fluid);
% Rewrite the transmissibilities
model.operators.T_all = T_all;
model.operators.T = T;

% The rock and model of Cartesian grid
rockC  = makeRock(GC, [.25, .25]*darcy, 0.25);
modelC = TwoPhaseOilWaterModel(GC, rockC, fluid);

%% Define wells for the hybrid grid
% Setup two injectors
% INJ1: south-west corner cell
% INJ2: north-east corner cell
D = sqrt( sum(G.cells.centroids.^2, 2) );
wcI = [find( D==min(D) ), find( D==max(D) )];
rate = 100*(meter^3/day);
W = [];
for j = 1 : 2
    W = addWell(W, G, rock, wcI(j), 'Name', sprintf('INJ%d',j), ...
        'sign', 1, 'comp_i', [1, 0], 'Val', rate, 'Type', 'rate');
end

% Producer in the hybrid grid
% Find the wellbore face first
D = bsxfun(@minus, G.faces.centroids, pW);
D = sqrt(sum(D.^2,2));
wf = find(~all(G.faces.neighbors,2) & D < 1.5*rW);
% The well index is the half transmissibility (radial) of the wellbore face
WI = T_all(wf);
% Get the well cells
wcP = sum(G.faces.neighbors(wf,:), 2);
bhp = 50*barsa;
W = addWell(W, G, rock, wcP, 'Name', 'PROD', 'sign', -1, ...
    'comp_i', [1, 1], 'Val', bhp, 'Type', 'bhp', 'WI', WI);

% Plot the well cells
figure, hold on; axis equal off
plotGrid(G, 'facecolor', 'none')
plotGrid(G, W(3).cells)
xlim([498,502]); ylim([498,502]);
legend('G', 'Well cells')

%% Define wells for the Cartesian grid
% Setup two injectors
% INJ1: south-west corner cell
% INJ2: north-east corner cell
D = sqrt( sum(GC.cells.centroids.^2, 2) );
wcI = [find( D==min(D) ), find( D==max(D) )];
WC = [];
for j = 1 : 2
    WC = addWell(WC, GC, rockC, wcI(j), 'Name', sprintf('INJ%d',j), 'sign', 1, ...
        'comp_i', [1, 0], 'Val', rate, 'Type', 'rate');
end

% Producer
D = bsxfun(@minus, GC.cells.centroids, pW);
D = sqrt(sum(D.^2,2));
wcP = find(D==min(D));
WC = addWell(WC, GC, rockC, wcP, 'Name', 'PROD', 'sign', -1, ...
    'comp_i', [1, 1], 'Val', bhp, 'Type', 'bhp', 'Radius', rW);

%% Define the schedule
timesteps = ones(20,1)*30*day;
schedule  = simpleSchedule(timesteps, 'W', W);
scheduleC = simpleSchedule(timesteps, 'W', WC);

%% Define the initial state
sW = zeros(G.cells.num, 1);
sat = [sW, 1 - sW];
state0 = initResSol(G, 300*barsa, sat);

sW = zeros(GC.cells.num, 1);
sat = [sW, 1 - sW];
state0C = initResSol(GC, 300*barsa, sat);

%% Run the simulation
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule);
[wellSolsC, statesC, reportC] = simulateScheduleAD(state0C, modelC, scheduleC);

%% Compare the well solutions
% Show the differences of producer's data (bhp, qWs, qOs) between the 
% hybrid grid and Cartesian grid
plotWellSols({wellSols, wellSolsC}, report.ReservoirTime)

%% Compare the pressure and oil saturation
ts = 20;

pargs = {'EdgeColor','none'};
figure
subplot(2,2,1), axis equal tight off
plotCellData(G, states{ts}.pressure/barsa, pargs{:})
title('Pressure of the hybrid grid')

subplot(2,2,2), axis equal tight off
plotCellData(GC, statesC{ts}.pressure/barsa, pargs{:})
title('Pressure of the Cartesian grid')

subplot(2,2,3), axis equal tight off
plotCellData(G, states{ts}.s(:,2), pargs{:})
title('Oil saturation of the hybrid grid')

subplot(2,2,4), axis equal tight off
plotCellData(GC, statesC{ts}.s(:,2), pargs{:})
title('Oil saturation of the Cartesian grid')

