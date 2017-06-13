% -------- Multi segment well example based on spe1-model ------------------

% Load necessary modules
mrstModule add ad-blackoil ad-props deckformat ad-core mrst-gui


%% Set up model based on spe1
[G, rock, fluid, deck, state] = setupSPE1();
[nx, ny] = deal(10);

model = selectModelFromDeck(G, rock, fluid, deck);

% Convert the deck schedule into a MRST schedule
schedule = convertDeckScheduleToMRST(model, deck);

%% Construct multisegment well, and replace producer in original example
% The well has 13 nodes 
% Nodes 1-7 represents the wellbore (frictional pressure drop)
% Well segments 2-8, 3-9, ..., 7-13 will represent valves
% Only nodes 8-13 are connected to reservoir (grid-cells c1-c6)

% 1 - 2 - 3 - 4 - 5 - 6 - 7
%     |   |   |   |   |   |
%     8   9  10  11  12  13
%     |   |   |   |   |   |
%    c1  c2  c3  c4  c5  c6

% connection grid cells (along edge of second layer)
c = nx*ny + (2:7)';
% First, initialize production well as "standard" well structure
prod = addWell([], G, rock, c, 'name', 'prod', 'refDepth', G.cells.centroids(1,3), ...
               'type', 'bhp', 'val', 250*barsa);

% Define additional properties needed for ms-well 
% We have 12 node-to-node segments
topo = [1 2 3 4 5 6 2 3  4  5  6  7
        2 3 4 5 6 7 8 9 10 11 12 13]';
%      |   tubing  |    valves     |     
% Create sparse cell-to-node mapping
cell2node = sparse((8:13)', (1:6)', 1, 13, 6);
% Segment lengths/diameters and node depths/volumes
lengths = [300*ones(6,1); nan(6,1)];
diam    = [.1*ones(6,1); nan(6,1)];
depths  = G.cells.centroids(c(1), 3)*ones(13,1);
vols    = ones(13,1);
% Convert to ms-well
prod = convert2MSWell(prod, 'cell2node', cell2node, 'topo', topo, 'G', G, 'vol', vols, ... 
                   'nodeDepth', depths, 'segLength', lengths, 'segDiam', diam);

% Finally, we set flow model for each segment:
%   segments  1-6: wellbore friction model
%   segments 7-12: nozzle valve model 
[wbix, vix]  = deal(1:6, 7:12);
roughness = 1e-4;
nozzleD   = .0025;
discharge = .7;
nValves   = 30; % number of valves per connection

prod.segments.flowModel = @(v, rho, mu)...
    [wellBoreFriction(v(wbix), rho(wbix), mu(wbix), prod.segments.diam(wbix), ...
                      prod.segments.length(wbix), roughness, 'massRate'); ...
     nozzleValve(v(vix)/nValves, rho(vix), nozzleD, discharge, 'massRate')];

% also define standard gas injector
inj = addWell([], G, rock, 100, 'name', 'inj', 'type', 'rate', 'Comp_i', [0 0 1], ...
              'val', 2.5e6/day,'refDepth', G.cells.centroids(1,3));

% combine inj and prod in struct W          
W = combineMSwithRegularWells(inj, prod);
schedule.control.W = W;
 
%% setup well-models and run simulation
wm = {SimpleWell(W(1)), MultisegmentWell(W(2))};

model.FacilityModel = model.FacilityModel.setupWells(W, wm);
model.extraStateOutput = true;

[wellSols, states, report] = simulateScheduleAD(state, model, schedule);

%% plotting
plotWellSols({wellSols}, report.ReservoirTime);
% plot pressure along wellbore for step 1 and step 120
figure, hold on
for k =  [1 120]
    plot([wellSols{k}(2).bhp; wellSols{k}(2).nodePressure(1:6)]/barsa, ...
         '-o', 'LineWidth', 2);
end
legend('Step 1', 'Step 120', 'Location', 'northwest')
set(gca, 'Fontsize', 14), xlabel('Well node'), ylabel('Pressure [bar]')

