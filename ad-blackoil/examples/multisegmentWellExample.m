%% Multi-segment well example based on SPE 1 benchmark model
% This example demonstrates the use of multisegment wells in MRST. The
% default wells in MRST assume flow in the well happens on a very short
% time-scale compared to the reservoir and can thus be accurately modelled
% by an instantaneous model.
%
% For some problems, however, more sophisticated well models may be
% required. The multisegment well class in MRST supports transient
% flow in the well. In addition, the class also allows for a more
% fine-grained representation of the well itself where the effect of valves
% and friction models can be included per segment.

% Load necessary modules
mrstModule add ad-blackoil ad-props deckformat ad-core mrst-gui example-suite

%% Set up model based on SPE 1
% We use the SPE 1 example to get a simple black-oil case to modify. For
% more details on this model, see the the example blackoilTutorialSPE1.m.
[G, rock, fluid, deck, state] = setupSPE1();
[nx, ny] = deal(G.cartDims(1),G.cartDims(2));

model = ThreePhaseBlackOilModel(G, rock, fluid, 'inputdata', deck);
% Add extra output
model.extraStateOutput = true;

% Convert the deck schedule into a MRST schedule
schedule = convertDeckScheduleToMRST(model, deck);

%% Construct multisegment well, and replace producer in original example
% The well has 13 nodes. Nodes correspond to volumes in the well and are
% analogous to the grid cells used when the reservoir itself is
% discretized. The well is perforated in six grid blocks and we add nodes
% both for the wellbore and the reservoir contacts. By modelling the nodes
% in the reservoir contacts and the wellbore separately, we can introduce
% a valve between the reservoir and the wellbore.
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

% Make a conceptual illustration of the multisegment horizontal well
clf
flag=true(G.cells.num,1); flag([2; c])=false;
plotGrid(G,flag,'FaceColor',[1 1 .7]);
plotGrid(G,c,'FaceColor','r','FaceAlpha',.3,'EdgeColor','r','LineWidth',1);
x = G.cells.centroids(c,[1 1 1]); x(:,2) = NaN;
y = G.cells.centroids(c,[2 2 2]); y(:,2) = NaN;
z = G.cells.centroids(c,[3 3 3]); z(:,2) = NaN;
z(:,1) = z(:,1)+1.5; z(:,3) = z(:,3)-1.5;
hold on
plot3(x(1,[1 1]), y(1,[1 1]), z(1,1) - [0 15], '-r', 'LineWidth',3);
plot3(x, y, z,'-or','MarkerSize',7,'MarkerFaceColor',[.6 .6 .6],'LineWidth',2);
plot3(x(:,[1 3 2])',y(:,[1 3 2])',z(:,[1 3 2])','-r','LineWidth',2);
hold off, view(-30,25), axis tight off

%%
% First, initialize production well as "standard" well structure
prodS = addWell([], G, rock, c, 'name', 'prod', ...
                'refDepth', G.cells.centroids(1,3), ...
               'type', 'rate', 'val', -8e5*meter^3/day);

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
depths  = G.cells.centroids(c([1 1:end 1:end]), 3);
vols    = ones(13,1);
% Convert to ms-well
prodMS = convert2MSWell(prodS, 'cell2node', cell2node, 'topo', topo, 'G', G, 'vol', vols, ...
                   'nodeDepth', depths, 'segLength', lengths, 'segDiam', diam);

% Finally, we set flow model for each segment:
%   segments  1-6: wellbore friction model
%   segments 7-12: nozzle valve model
% The valves have 30 openings per connection
[wbix, vix]  = deal(1:6, 7:12);
[roughness, nozzleD, discharge, nValves] = deal(1e-4, .0025, .7, 30);

% Set up flow model as a function of velocity, density and viscosity
prodMS.segments.flowModel = @(v, rho, mu)...
    [wellBoreFriction(v(wbix), rho(wbix), mu(wbix), prodMS.segments.diam(wbix), ...
                      prodMS.segments.length(wbix), roughness, 'massRate'); ...
     nozzleValve(v(vix)/nValves, rho(vix), nozzleD, discharge, 'massRate')];

% In addition, we define a standard gas injector. Different well types can
% easily be combined in MRST, each with their own models for pressure drop
% and solution variables.
inj = addWell([], G, rock, 100, 'name', 'inj', 'type', 'rate', 'Comp_i', [0 0 1], ...
              'val', 2.5e6*meter^3/day,'refDepth', G.cells.centroids(1,3));

%% Run schedule with simple wells
% We first simulate a baseline where the producer is treated as a simple
% well with instantaneous flow and one node per well-reservoir contact.
W_simple = [inj; prodS];
schedule.control.W = W_simple;

[wellSolsSimple, statesSimple] = simulateScheduleAD(state, model, schedule);

%% Set up well models
% We combine the simple and multisegment well together
W = combineMSwithRegularWells(inj, prodMS);
schedule.control.W = W;
% Initialize the facility model. This is normally done automatically by the
% simulator, but we do it explicitly on the outside to view the output
% classes. This approach is practical if per-well adjustments to e.g.
% tolerances are desired.

% First, validate the model to set up a FacilityModel
model = model.validateModel();
% Then apply the new wells to the FacilityModel
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
figure, hold all
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
            'datasetnames', {'Multisegment', 'Standard'});
        
%% Show the advancing displacement front
mrstModule add coarsegrid
pg = generateCoarseGrid(G,ones(G.cells.num,1));
figure, h = [];
plotGrid(pg,'FaceColor','none');
plotGrid(G,c,'FaceColor','none','EdgeColor','r');
plotWell(G,[inj; prodS]);
 caxis([0 .5]); view(37.5,34); caxis([0 .5]); axis tight; colorbar, drawnow;
for i=1:numel(states)
    sg = states{i}.s(:,3); inx = sg>1e-5;
    if sum(inx)>0
        delete(h)
        h=plotCellData(G,sg,sg>1e-4,'EdgeAlpha',.01);
        colorbar
    drawnow;
    end
end

%% Visualize drawdown in well as function of time
% We extract the pressure in all the nodes and visualize this, together
% with the pressure in the completed reservoir cells as function of time
[xx,tt]=meshgrid(1:6,report.ReservoirTime/day);
pr = cellfun(@(x) x.pressure(c)', states, 'UniformOutput',false);
pa = cellfun(@(x) x(2).nodePressure(7:12)', wellSols, 'UniformOutput',false);
pw = cellfun(@(x) x(2).nodePressure(1:6)', wellSols, 'UniformOutput',false);
figure
hold on
surfWithOutline(xx, tt, vertcat(pr{:})/barsa);
surfWithOutline(xx, tt, vertcat(pa{:})/barsa);
surfWithOutline(xx, tt, vertcat(pw{:})/barsa);
hold off, box on, axis tight
view(50,10), shading interp
set(gca,'Projection','Perspective');
akse = axis();
cax = caxis();

%%
pr = cellfun(@(x) x.pressure(c)', statesSimple, 'UniformOutput',false);
pw = cellfun(@(x) x(2).bhp + x(2).cdp', wellSolsSimple, 'UniformOutput',false);
figure
hold on
surfWithOutline(xx, tt, vertcat(pr{:})/barsa);
surfWithOutline(xx, tt, vertcat(pw{:})/barsa);
hold off, box on, axis tight, axis(akse), caxis(cax)
view(50,10), shading interp
set(gca,'Projection','Perspective');

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
