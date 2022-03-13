%% Simulation of an inverted five-spot pattern with dG
% In this example, we simulate water injection in an inverted five-spot
% pattern posed on a PEBI mesh with dG(0) and dG(1) and visualize the
% results

%% Add modules
mrstModule add dg upr vem vemmech spe10 ad-props ad-core ad-blackoil ...
    sequential

%% Set up fine-scale model
% We extract the lower half of layer 13 of SPE10 2
[state0Ref, modelRef, scheduleRef] = setupSPE10_AD('layers', 13           , ...
                                                   'J'     , (1:110) + 30);
GRef    = modelRef.G;
rockRef = modelRef.rock;
WRef    = scheduleRef.control(1).W;
fluid = modelRef.fluid;

%% Make PEBI grid
% We use the upr module to construct a PEBI grid with refinement around the
% wells
rng(2019) % For reproducibility
n = 20;   % Approx number of cells in each direction
% Get well coordinates
l = max(GRef.nodes.coords(:,1:2));
wellLines = mat2cell(GRef.cells.centroids(vertcat(WRef.cells),1:2), ...
                                                    ones(numel(WRef),1), 2)';
% Construct PEBI grid
G = pebiGrid2D(max(l)/n, l, 'cellConstraints', wellLines, ... % Well coords
                            'CCRefinement'   , true     , ... % Refine
                            'CCFactor'       , 0.4      );
G = computeGeometry(G);       % Compute geometry
G = computeCellDimensions(G); % Compute cell dimensions

%% Sample rock properties
% We assign rock properties in the PEBI grid cells by sampling from the
% fine grid using sampleFromBox
poro = sampleFromBox(G, reshape(rockRef.poro, GRef.cartDims));
perm = zeros(G.cells.num,G.griddim);
for i = 1:G.griddim
    perm(:,i) = sampleFromBox(G, reshape(rockRef.perm(:,i), GRef.cartDims));
end
rock = makeRock(G, perm, poro);

%% Set up schedule and initial state
W        = WRef;
schedule = scheduleRef;
x  = G.cells.centroids;
xwR = GRef.cells.centroids(vertcat(WRef.cells),1:2);
% Slick oneliner to find corresponding cells in the new grid
[~, c] = min(sum(bsxfun(@minus, reshape(xwR, [], 1 , G.griddim), ...
                             reshape(x , 1 , [], G.griddim)).^2,3), [], 2);
for i = 1:numel(W)
    W(i).cells = c(i);
end
xw = G.cells.centroids(vertcat(W.cells),1:2);
schedule.control.W = W;
state0 = initResSol(G, state0Ref.pressure(1), state0Ref.s(1,:));

%% Inspect the geological model
figure('Position', [0,0,800,410]), subplot(1,2,1)
Kref = convertTo(rockRef.perm(:,1),milli*darcy);
plotCellData(GRef,log10(Kref),'edgeAlpha',.1); axis tight
set(gca,'FontSize',12)
mrstColorbar(Kref,'South',true); cx = caxis();
hold on; plot(xwR(:,1),xwR(:,2),'.r','MarkerSize',18); hold off

subplot(1,2,2);
K = convertTo(rock.perm(:,1),milli*darcy);
plotCellData(G,log10(K),'edgeAlpha',.1); axis tight; caxis(cx);
set(gca,'FontSize',12)
mrstColorbar(K,'South',true); axis tight
hold on; plot(xw(:,1),xw(:,2),'.r','MarkerSize',18); hold off

%% Set base model
% The base model is a two-phase oil-water model
model        = GenericBlackOilModel(G, rock, fluid, 'gas', false);
pmodel       = PressureModel(model); % Pressure model
makeSeqModel = @(tmodel) ...
    SequentialPressureTransportModel(pmodel, tmodel,'parentModel', model);

%% Simulate with dG(0)
tmodelDG0 = TransportModelDG(model, 'degree', 0);
modelDG0  = makeSeqModel(tmodelDG0);
[wsDG0, stDG0, repDG0] = simulateScheduleAD(state0, modelDG0, schedule);

%% Simulate with dG(1)
tmodelDG1 = TransportModelDG(model, 'degree', 1);
modelDG1  = makeSeqModel(tmodelDG1);
[wsDG1, stDG1, repDG1] = simulateScheduleAD(state0, modelDG1, schedule);

%% Visualize results
% We plot the evolving saturation front computed using the two
% discuretizations as surface plots. To this end, we use plotSaturationDG,
% which supports surface plots of cell-wise discontinuous data using the
% patch function.
% Axis properties
setAxProps = @(ax) set(ax, 'Projection'        , 'Perspective', ...
                           'PlotBoxAspectRatio', [1,1,0.3]    , ...
                           'View'              , [-75,52]     , ...
                           'XLim'              , [0,l(1)]     , ...
                           'YLim'              , [0,l(2)]     , ...
                           'ZLim'              , [0,1]        , ...
                           'Box'               , 'on'         );
% For plotting wells
pw = @() plotWell(GRef, WRef, 'height', -1, 'color', 'k');
% Get coordinates for plotting (to be used in patch)
coords = getPlotCoordinates(G);
close all; figure('Position', [0,0,1000,500]);
subplot(1,2,1); % Plot dG(0)
[hs0, satDG0] = plotSaturationDG(tmodelDG0.discretization, stDG0{1}, ...
                                                         'coords', coords);
pw(); setAxProps(gca); caxis([0.2,0.8]); camlight; % Set axis properties
subplot(1,2,2); % Plot dG(1)
[h1, satDG1] = plotSaturationDG(tmodelDG1.discretization, stDG1{1}, ...
                                                         'coords', coords);
pw(); setAxProps(gca); caxis([0.2,0.8]); camlight; % Set axis properties

for i = 1:numel(stDG1)
    % Update dG(0) patch
    s0 = satDG0(stDG0{i});
    hs0.Vertices(:,3)   = s0;
    hs0.FaceVertexCData = s0;
    % Update dG(1) patch
    s1 = satDG1(stDG1{i});
    h1.Vertices(:,3)   = s1;
    h1.FaceVertexCData = s1;
    % Pause
    pause(0.2),
end

%% Plot the difference between the two solutions
figure
plotCellData(G,stDG1{end}.s(:,1)-stDG0{end}.s(:,1),'EdgeColor','none'); axis tight
view(-90,90), box on
colormap(interp1([0; 0.5; 1], [1, 0, 0; 1, 1, 1; 0, 0, 1], 0:0.01:1))
caxis([-.151 .151]);
colorbar('EastOutside');

%% Plot the production responses
plotWellSols({wsDG0,wsDG1},'datasetnames',{'dG(0)','dG(1)'},'field','qWs');

%% Plot the number of iterations taken in different parts of the solver
reports = {repDG0, repDG1};
ns = numel(reports);
stats = cell(1, ns);
for i = 1:ns
    stats{i} = getPressureTransportIterations(reports{i});
end
total = zeros(ns, 3);
for i = 1:ns
    s = stats{i};
    total(i, 1) = sum(s.pressure);
    total(i, 2) = sum(s.transport);
    total(i, 3) = sum(s.outer);
end
figure;
bar(total)
set(gca, 'XTickLabel', {'dG(0)','dG(1)'}, 'XTickLabelRotation',20);
legend('Pressure', 'Transport', 'Outer','Location','best');
text(1:ns,total(:,2),num2str(total(:,2)),'vert','bottom','horiz','center');

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
