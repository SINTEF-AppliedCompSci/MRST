%% Gas injection problem
% In this example, we consider the same setup as in the two-phase example
% with wells, except that we now assume a three-phase, compressible model
% with gas injection. For comparison, we also simulate injection in the
% same reservoir geometry without (embedded) fractures. The purpose of the
% example is to demonstrate how the methods from the HFM module easily can
% be combined with solvers from the ad-blackoil module.

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add ad-props ad-core % AD framework
mrstModule add ad-blackoil      % Three phase simulator
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

%% Grid and fracture lines
% Construct a Cartesian grid comprising 50-by-20 cells, where each cell has
% dimension 10-by-10 m^2. Define 3 fracture lines by supplying their end
% points in the form [x1 y1 x2 y2].

celldim = [50 20];
physdim = [500 200];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [80,  100, 270, 180;...
      130, 160, 340,  40;...
      260,  40, 420, 150] ; % fractures lines in [x1 y1 x2 y2] format.

%% Process fracture lines
% Using the input fracture lines, we identify independent fracture networks
% comprising of connected lines. In this example, there is only 1 fracture
% network consisting of 3 fracture lines. We also identify the fine-cells
% in the matrix containing these fractures. Fracture aperture is set to
% 0.04 meters. The matrix grid and fracture lines are plotted.

dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracture] = processFracture2D(G, fl, 'verbose', mrstVerbose);
fracture.aperture = 1/25; % Fracture aperture
figure;
plotFractureLines(G,fracture);
axis equal tight; 
box on

%% Compute CI and construct fracture grid
% For each matrix block containing a fracture, we compute a fracture-matrix
% conductivity index (CI) in accordance with the hierarchical fracture
% model (a.k.a. embedded discrete fracture model). Following that, we
% compute a fracture grid where the fracture cell size is defined to be
% 10 m. Next, the fracture grid is plotted on top of the matrix grid.

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 5; cell_size = 10; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); 
axis equal tight; box on

%% Set rock properties in fracture and matrix
% For this example, we will generate the porosity as a Gaussian field. To
% get a crude approximation to the permeability-porosity relationship, we
% assume that our medium is made up of uniform spherical grains of diameter
% dp = 10 m, for which the specic surface area is Av = 6 = dp. With these
% assumptions, we can then use the Carman Kozeny relation to calculate the
% isotropic permeability (K). The rock properties are then plotted. We also
% identify fracture-matrix and fracture-fracture connections and compute
% transmissibility for each connection.

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
p = gaussianField(celldim, [0.2 0.4], [11 3], 2.5);
mp = 0.01;
p(p < mp) = mp;
K = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);

p = p(:);
K = K(:);

G.rock.poro = p(:);
G.rock.perm = K(:);
K_frac = 10000; % Darcy
G = makeRockFrac(G, K_frac, 'porosity', 0.8);

clf; plotToolbar(G, G.rock); axis tight equal
line(fl(:,1:2:3)',fl(:,2:2:4)',1e-3*ones(2,size(fl,1)),'Color','r','LineWidth',0.5);

[G,T] = defineNNCandTrans(G,F,fracture);


%% Define fluid properties
% Define a two-phase fluid model without capillarity. The fluid model has
% the values:
%
% * densities: [rho_w, rho_o] = [1000 700 250] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 5 0.2] cP.
% * corey-coefficient: [2, 2] = [2 2 2].

fluid = initSimpleADIFluid('mu' , [   1,  5, 0.2] .* centi*poise     , ...
                           'rho', [1000, 700, 250] .* kilogram/meter^3, ...
                           'n'  , [   2,   2, 2]);

pRef = 100*barsa;
c_w = 1e-8/barsa;
c_o = 1e-5/barsa;
c_g = 1e-3/barsa;

fluid.bW = @(p) exp((p - pRef)*c_w);
fluid.bO = @(p) exp((p - pRef)*c_o);
fluid.bG = @(p) exp((p - pRef)*c_g);

%% Define three-phase compressible flow model
% We define a three-phase black-oil model without dissolved gas or vaporized
% oil. This is done by first instantiating the blackoil model, and then
% manually passing in the internal transmissibilities and the topological
% neighborship from the embedded fracture grid. This will cause a warning
% that can simply be ignored.

model = ThreePhaseBlackOilModel(G, [], fluid, 'disgas', false, 'vapoil', false);
N = getNeighbourship(G, 'topological', true);
intx = all(N ~= 0, 2);

% Send in internal transmissibility and neighborship to the operator setp

model.operators = setupOperatorsTPFA(G, G.rock, 'trans', T(intx), 'neighbors', N(intx, :));

%% Add wells
% We have a single gas injector that injects a total of one pore volume of
% gas at surface conditions over a period of five years. Two producers are
% added and set to bottom hole pressure controls of 50 bar.

totTime = 5*year;
tpv = sum(model.operators.pv);

inj = 1;
prod1 = celldim(1)*celldim(2);
prod2 = celldim(1);
wellRadius = 0.1;

% Injector
W = addWell([], G.Matrix, G.Matrix.rock, inj,'InnerProduct', 'ip_tpf','Type', ...
    'rate', 'Val', tpv/totTime, 'Radius', wellRadius, 'Comp_i', [0, 0, 1]);
% First producer
W = addWell(W, G.Matrix, G.Matrix.rock, prod1, 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 50*barsa, 'Radius', wellRadius, 'Comp_i', [1, 1, 1]/3);
% Second producer
W = addWell(W, G.Matrix, G.Matrix.rock, prod2, 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 50*barsa, 'Radius', wellRadius, 'Comp_i', [1, 1, 1]/3);

%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water and oil initially. We also set up a simple-time step strategy that
% ramps up gradually towards 30 day time-steps.

s0 = [0.2, 0.8, 0];
state  = initResSol(G, pRef, s0);
dt = rampupTimesteps(totTime, 30*day, 8);
schedule = simpleSchedule(dt, 'W', W);

%% Simulate problem
%fn = getPlotAfterStep(state, model, schedule);
[ws, states, report] = simulateScheduleAD(state, model, schedule, ...
   'afterStepFn', getPlotAfterStep(state, model, schedule));

%% Create and solve the same problem without any fractures present
% Fractures can have a large impact on fluid displacement, but are often
% associated with great uncertainty in terms of locations, orientation and
% permeability. To demonstrate the effect of fractures on the transport, we
% create another problem where the fractures themselves have been omitted.

G0 = cartGrid(celldim, physdim);
G0 = computeGeometry(G0);
% Create rock without fractures
rock0 = makeRock(G0, K, p);
% Same model, but different grid and rock
model0 = ThreePhaseBlackOilModel(G0, rock0, fluid, 'disgas', false, 'vapoil', false);
state0  = initResSol(G0, pRef, s0);

% Make a copy of the wells in the new grid
W0 = [];
for i = 1:numel(W)
    w = W(i);
    W0 = addWell(W0, G0, rock0, w.cells, 'type', w.type, 'val', w.val, ...
                 'wi', w.WI, 'comp_i', w.compi);
end
% Set up and simulate the schedule
schedule0 = simpleSchedule(dt, 'W', W0);
[ws0, states0, report0] = simulateScheduleAD(state0, model0, schedule0, ...
   'afterStepFn', getPlotAfterStep(state0, model0, schedule0));

%% Plot the results
% We plot the production curves for the two wells and compare between the
% two cases, as well as the final gas saturation.

tm = report.ReservoirTime/day;
flds = {'qGs', 'qWs', 'qOs'};
names = {'Gas production', 'Water Production', 'Oil production'};
sgn = vertcat(W.sign);
inj = find(sgn > 0);
prod = find(sgn <= 0);
colors = lines(nnz(prod));

for i = 1:numel(flds)
    figure(i); clf
    q0 = -getWellOutput(ws0, flds{i}, prod)*day;
    q =  -getWellOutput(ws, flds{i}, prod)*day;
    l = {};
    for j = 1:numel(prod)
        c = colors(j, :);
        
        plot(tm, q(:, j), 'color', c);
        hold on
        plot(tm, q0(:, j), '--', 'color', c, 'linewidth', 2);
        title(names{i})
        ylabel('Production [m^3/day]');
        xlabel('Time [days]')
        n = W(prod(j)).name;
        l = [l, {[n, ' with fracture'], [n, ' without fracture']}]; %#ok
    end
    legend(l)
end

figure;
subplot(2, 1, 1)
plotCellData(G0, states{end}.s(1:G0.cells.num, 3), 'EdgeColor', 'none')
line(fl(:,1:2:3)',fl(:,2:2:4)',1e-3*ones(2,size(fl,1)),'Color','r','LineWidth',0.5);
axis equal tight
caxis([0, 1])
title('Gas saturation with fractures')
subplot(2, 1, 2)
plotCellData(G0, states0{end}.s(1:G0.cells.num, 3), 'EdgeColor', 'none')
title('Gas saturation without fractures')
axis equal tight
caxis([0, 1])
xw = G.cells.centroids(vertcat(W.cells),:);
for i=1:2,
   subplot(2,1,i)
   hold on
   plot3(xw(1,1),xw(1,2),1e-3,'xk','LineWidth',2,'MarkerSize',8);
   plot3(xw(2,1),xw(2,2),1e-3,'ok','LineWidth',2,'MarkerSize',8);
   plot3(xw(3,1),xw(3,2),1e-3,'ok','LineWidth',2,'MarkerSize',8);
   hold off
end

%% Interactive plotting
% Compare the well cruves

plotWellSols({ws, ws0}, report.ReservoirTime, 'datasetnames', {'With fractures', 'Without fractures'})

figure;
plotToolbar(G, states)
line(fl(:,1:2:3)',fl(:,2:2:4)',1e-3*ones(2,size(fl,1)),'Color','r','LineWidth',0.5);
hold on
plot3(xw(1,1),xw(1,2),1e-3,'xk','LineWidth',2,'MarkerSize',8);
plot3(xw(2,1),xw(2,2),1e-3,'ok','LineWidth',2,'MarkerSize',8);
plot3(xw(3,1),xw(3,2),1e-3,'ok','LineWidth',2,'MarkerSize',8);
hold off
axis equal tight
title('With fractures')

figure;
plotToolbar(G0, states0)
hold on
plot3(xw(1,1),xw(1,2),1e-3,'xk','LineWidth',2,'MarkerSize',8);
plot3(xw(2,1),xw(2,2),1e-3,'ok','LineWidth',2,'MarkerSize',8);
plot3(xw(3,1),xw(3,2),1e-3,'ok','LineWidth',2,'MarkerSize',8);
hold off
axis equal tight
title('Without fractures')

% <html>
% <p><font size="-1">
% Copyright 2009-2016 TU Delft and SINTEF ICT, Applied Mathematics.
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