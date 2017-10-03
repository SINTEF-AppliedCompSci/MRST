%% Multi-segment well example based on SPE 1 benchmark model
% This example demonstrates the use of multisegment wells in MRST. The
% default wells in MRST assume flow in the well happens on a very short
% time-scale compared to the reservoir and can thus be accurately modelled
% by a instantaneous model. 
%
% For some problems, however, more sophisticated well models may be
% required. The multisegment well class in MRST supports transient
% flow in the well. In addition, the class also allows for a more
% fine-grained representation of the well itself where the effect of valves
% and friction models can be included per segment.

% Load necessary modules
mrstModule add ad-blackoil ad-props deckformat ad-core mrst-gui

%% Set up model based on SPE 1
% We use the SPE 1 example to get a simple black-oil case to modify. For
% more details on this model, see the the example blackoilTutorialSPE1.m.
[G, rock, fluid, deck, state] = setupSPE1();
[nx, ny] = deal(10);

model = selectModelFromDeck(G, rock, fluid, deck);
% Add extra output 
model.extraStateOutput = true;

% Convert the deck schedule into a MRST schedule
schedule = convertDeckScheduleToMRST(model, deck);

%% Construct multisegment well, and replace producer in original example
% The well has 13 nodes. Nodes correspond to volumes in the well and are
% analogous to the grid cells used when the reservoir itself is
% discretized. The well is perforated in six grid blocks and we add nodes
% both for the well bore and the reservoir contacts. By modelling the nodes
% in the reservoir contacts and the well bore separately, we can introduce
% a valve between the reservoir and the well bore. 
%
% Nodes 1-7 represents the wellbore (frictional pressure drop)
% Well segments 2-8, 3-9, ..., 7-13 will represent valves
% Only nodes 8-13 are connected to reservoir (grid-cells c1-c6)
%
%   Segments with friction model for pressure drop
%   |   |   |   |   |   |
%   v   v   v   v   v   v
% 1 - 2 - 3 - 4 - 5 - 6 - 7  <- Nodes in well bore
%     |   |   |   |   |   |  <- Segments with valves for pressure drop
%     8   9  10  11  12  13  <- Nodes between reservoir and wellbore
%     |   |   |   |   |   |  <- Pressure drop by standard well index model
%    c1  c2  c3  c4  c5  c6  <- Reservoir cells

% connection grid cells (along edge of second layer)
c = nx*ny + (2:7)';
% First, initialize production well as "standard" well structure
prod0 = addWell([], G, rock, c, 'name', 'prod', 'refDepth', G.cells.centroids(1,3), ...
               'type', 'bhp', 'val', 250*barsa);
prod = prod0;
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
% Set up flow model as a function of velocity, density and viscosity
prod.segments.flowModel = @(v, rho, mu)...
    [wellBoreFriction(v(wbix), rho(wbix), mu(wbix), prod.segments.diam(wbix), ...
                      prod.segments.length(wbix), roughness, 'massRate'); ...
     nozzleValve(v(vix)/nValves, rho(vix), nozzleD, discharge, 'massRate')];

% In addition, we define a standard gas injector. Different well types can
% easily be combined in MRST, each with their own models for pressure drop
% and solution variables.
inj = addWell([], G, rock, 100, 'name', 'inj', 'type', 'rate', 'Comp_i', [0 0 1], ...
              'val', 2.5e6/day,'refDepth', G.cells.centroids(1,3));

%% Run schedule with simple wells
% We first simulate a baseline where the producer is treated as a simple
% well with instantaneous flow and one node per well-reservoir contact.
W_simple = [inj; prod0];
schedule.control.W = W_simple;

[wellSolsSimple, statesSimple] = simulateScheduleAD(state, model, schedule);

%% Set up well models
% We combine the simple and multisegment well together
W = combineMSwithRegularWells(inj, prod);
schedule.control.W = W;
% Initialize the facility model. This is normally done automatically by the
% simulator, but we do it explicitly on the outside to view the output
% classes. This approach is practical if per-well adjustments to e.g.
% tolerances are desired.
model.FacilityModel = model.FacilityModel.setupWells(W);

% View the simple injector
disp(model.FacilityModel.WellModels{1})
% View the multisegment producer
disp(model.FacilityModel.WellModels{2})

%% Run the same simulation with multisegment wells
% Note that if verbose output is enabled, the convergence reports will
% include additional fields corresponding to the well itself.
[wellSols, states, report] = simulateScheduleAD(state, model, schedule);
%% Plot the well-bore pressure in the multisegment well
% Plot pressure along wellbore for step 1 and step 120 (final step)
figure, hold on
for k =  [1 120]
    plot([wellSols{k}(2).bhp; wellSols{k}(2).nodePressure(1:6)]/barsa, ...
         '-o', 'LineWidth', 2);
end
legend('Step 1', 'Step 120', 'Location', 'northwest')
set(gca, 'Fontsize', 14), xlabel('Well node'), ylabel('Pressure [bar]')

%% Launch interactive plotting
% We can compare the two different well models interactively and observe
% the difference between the two modelling approaches.
plotWellSols({wellSols, wellSolsSimple}, report.ReservoirTime, ...
            'datasetnames', {'Complex wells', 'Standard wells'});
