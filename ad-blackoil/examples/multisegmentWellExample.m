% -------- Multi segment well example based on spe1-model ------------------
% Well topology hard-coded and "standard" well converted using convert2MSWell.m

% Load necessary modules
mrstModule add ad-blackoil ad-props deckformat ad-core mrst-gui


%% Set up model based on spe1
[G, rock, fluid, deck, state] = setupSPE1();
[nx, ny] = deal(10);

model = selectModelFromDeck(G, rock, fluid, deck);

% Convert the deck schedule into a MRST schedule
schedule = convertDeckScheduleToMRST(model, deck);

%% Construct multisegment well from scratch, and replace producer in original example
% The well has 14 nodes as shown below (the top-node has index 0)
% Nodes 1-7 represents the wellbore (frictional pressure drop)
% Well segments 8-2, 9-3, ..., 13-7 will represent valves
% Only nodes 8-13 are connected to reservoir (gridcells c1-c6)

% 0
% |
% 1 - 2 - 3 - 4 - 5 - 6 - 7
%     |   |   |   |   |   |
%     8   9  10  11  12  13
%     |   |   |   |   |   |
%    c1  c2  c3  c4  c5  c6

% connection grid cells (along edge of second layer)
c = nx*ny + (2:7)';
% combine in initially in a "standard" well structure
prod = addWell([], G, rock, c, 'name', 'prod', 'refDepth', G.cells.centroids(1,3), ...
               'type', 'bhp', 'val', 250*barsa);

% well topology (node-to-node segments), ommit special top segments
topo = [1 2 3 4 5 6 2 3  4  5  6  7
        2 3 4 5 6 7 8 9 10 11 12 13]';
%      |   tubing  |    valves     |     
% create sparse cell-to-node mapping
cell2node = sparse((8:13)', (1:6)', 1, 13, 6);
% segment lengths/diameters and node depths/volumes
lengths = [300*ones(6,1); nan(6,1)];
diam    = [.1*ones(6,1); nan(6,1)];
depths  = G.cells.centroids(c(1), 3)*ones(13,1);
vols    = ones(13,1);

prod = convert2MSWell(prod, 'cell2node', cell2node, 'topo', topo, 'G', G, 'vol', vols, ... 
                   'nodeDepth', depths, 'segLength', lengths, 'segDiam', diam);

% set flow model for each segment. Friction for first six, valve for last six
[wbix, vix]  = deal(1:6, 7:12);
roughness = 1e-4;
nozzleD   = .0025;
discharge = .7;
nValves   = 30; % number of valves per connection

prod.segments.flowModel = @(v, rho, mu)...
    [wellBoreFriction(v(wbix), rho(wbix), mu(wbix), prod.segments.diam(wbix), ...
                      prod.segments.length(wbix), roughness, 'massRate'); ...
     ICDFriction(v(vix)/nValves, rho(vix), nozzleD, discharge, 'massRate')];

% also define standard gas injector
inj = addWell([], G, rock, 100, 'name', 'inj', 'type', 'rate', 'Comp_i', [0 0 1], 'val', 2e6/day,'refDepth', G.cells.centroids(1,3));
 

 W = combineMSwithRegularWells(inj, prod);
schedule.control.W = W;
 

wm = {SimpleWell(W(1)), MultisegmentWell(W(2))};

model.FacilityModel = model.FacilityModel.setupWells(W, wm);
model.extraStateOutput = true;


nSteps = 120;
schedule.step.val       = schedule.step.val(1:min(nSteps, 120));
schedule.step.control   = schedule.step.control(1:min(nSteps, 120));

[wellSols, states, report] = simulateScheduleAD(state, model, schedule);

