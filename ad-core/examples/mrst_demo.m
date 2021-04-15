%% MRST demo
% This example will demonstrate how to build and run a two-phase black-oil 
% simulation from scratch in MRST. This example is equivalent to the
% mrst_demo.m script shown in  the talk by Francesca Watson at the MATLAB 
% Energy Conference 2020: "MATLAB Reservoir Simulation Toolbox in Action"
% https://se.mathworks.com/videos/matlab-reservoir-simulation-toolbox-in-action-1610451936890.html
% 
% To run the example (and most other examples in MRST) it is best to run in
% cell mode and advance through the script step-by-step. 
%
% Inspecting individual functions will give more information on the input
% parameters to the function in question and the resulting outputs.


%% Add modules
% First we add the required mrst modules which is effectively the same as
% adding the module folders to the matlab path.
% - ad stands for automatic differentiation
% - ad-core and ad-blackoil add the simulator framework and the blackoil
% module. 
% - ad-props contains functionality for defining fluid properties.
% - mrst-gui contains functionality for interactive visualisation of 
% results.
mrstModule add ad-core ad-blackoil ad-props mrst-gui

% List the loaded modules.
mrstModule list

%% Build horizon for top structure
% Here we build our own grid by creating horizons and interpolating between
% horizons to get the corner-point grid in mrst format.

% Define areal mesh
[xmax,ymax, n]  = deal(1000*meter, 1000*meter, 30);
[x, y] = meshgrid(linspace(0,xmax,n+1), linspace(0,ymax,n+1));
[x, y] = deal(x',y');

% Basic dome structure
dome = 1-exp(sqrt((x - xmax/2).^2 + (y - ymax/2).^2)*1e-3);

% Periodic and random small-scale perturbation
[xn,yn] = deal(pi*x/xmax,pi*y/ymax);
perturb = sin(5*xn) + .5*sin(4*xn+6*yn) + cos(.25*xn*yn./pi^2) + cos(3*yn);
perturb = perturb/3.5;
rng(0);
[h, hr] = deal(8,1);
zt = 50 + h*perturb + rand(size(x))*hr - 20*dome;

% Visualise
surf(x,y,zt-.2, 'EdgeC','r','FaceC',[.8 .8 .8]),  hold on
set(gca,'ZDir','reverse')
view(-50,10); axis off

%% Build remaining horizons and create horizon structure
zb = zt + 30;
zmb = min(zb + 4 + 0.01*x - 0.020*y + hr*rand(size(x)), zb);
zmt = max(zb -15 + 0.01*x - 0.025*y + hr*rand(size(x)), zt);

horizons = {struct('x', x, 'y', y, 'z', zt), ...
    struct('x', x, 'y', y, 'z', zmb), ...
    struct('x', x, 'y', y, 'z', zb)};

% Visualise
clf
ht = surf(x,y,zt, 'EdgeC',[.5 .5 .5],'FaceC',[1 .4 .4]);  hold on;
hm = surf(x,y,zmb,'EdgeC',[.5 .5 .5],'FaceC',[.4 .4 1]); 
hb = surf(x,y,zb, 'EdgeC',[.5 .5 .5],'FaceC',[.4 1 .4]); hold off
set(gca,'ZDir','reverse')
view(-50,10); axis off; camlight


%% Interpolate to build unfaulted corner-point grid
% We first specify the horizonal dimensions of the grid and the number of
% layers between each horizon.
dims   = [40, 40]; 
layers = [6 3];
% Then we create the corner-point grid and finally convert to an mrst grid
% structure, called G by covention.
grdecl = convertHorizonsToGrid(horizons, 'dims', dims, 'layers', layers);
G      = processGRDECL(grdecl);

% Visualise
hold on
set([ht,hm,hb],'EdgeColor','none');
plotGrid(G,'FaceColor','none','FaceAlpha',.2,'EdgeAlpha',.5);


%% Insert faults
% We can create faults in the grid by modifiying the coordinates of the
% corner-points of certain cells and recreating the mrst grid.
[X,Y,Z]  = buildCornerPtNodes(grdecl);

i=47:80; Z(i,:,:) = Z(i,:,:) + .022*min(0,Y(i,:,:)-550);
j= 1:30; Z(:,j,:) = Z(:,j,:) + .021*min(0,X(:,j,:)-400);
j=57:80; Z(:,j,:) = Z(:,j,:) + .023*min(0,X(:,j,:)-750);
grdecl.ZCORN = Z(:);

G = processGRDECL(grdecl);
G = computeGeometry(G);
[~,~,k] = gridLogicalIndices(G);

% Visualise
clf, plotCellData(G,k,'EdgeAlpha',.2); view(3);
plotFaces(G,find(G.faces.tag>0),'EdgeColor','r','FaceColor',[.8 .8 .8]);
colormap(.7*jet + .3*ones(size(jet)));
view(-50,30); axis tight off

%% Petrophysics
% Set up permeability based on K-indices and introduce anisotropy by
% setting K_z = .1*K_x
% We use 4 layers with log-normal permeability distributions which have
% means values of [100 400 10 50]*milli*darcy respectively.

rng(357371);
[K,L] = logNormLayers(G.cartDims, [100 400 10 50]*milli*darcy);
K = K(G.cells.indexMap);
perm = [K, K, 0.1*K];
rock = makeRock(G, perm, 0.3);

% Visualise
clf;
K = convertTo(K,milli*darcy);
plotCellData(G, log10(K),'EdgeAlpha',.1);
mrstColorbar(K,'South',true,[1 1500]);
view(-50, 50), axis tight off

%% Define wells
% The model contains three production wells, operating at fixed bottom-hole
% pressure and one injection well with fixed injection rate, perforated 
% throughout all layers of the model. 

% Producers
simTime = 10*year;
pv      = poreVolume(G, rock);
injRate = 1*sum(pv)/simTime;
offset  = 10;

W = verticalWell([], G, rock, offset, offset, [],...
                'Name', 'P1', 'comp_i', [1 0], ...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50);
W = verticalWell(W, G, rock,  offset, floor(G.cartDims(1)/2)+3, [],...
                'Name', 'P2', 'comp_i', [1 0], ...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50);
W = verticalWell(W, G, rock, offset, G.cartDims(2) - offset/2, [], ...
                'Name', 'P3', 'comp_i', [1 0], ...
                'Val', 250*barsa, 'Type', 'bhp', 'refDepth', 50);

% Injectors
W = verticalWell(W, G, rock, G.cartDims(1)-5, offset, [],...
                'Name', 'I1', 'comp_i', [1 0], ...
                'Val', injRate, 'Type', 'rate', 'refDepth', 50);

% Visualise
plotWell(G, W,'color','k')
axis tight


%% Two-phase fluid model
% We setup a simple two phase fluid model with water and oil. For more
% information about the input parameters have a look at
% initSimpleADIFluid.m
fluid = initSimpleADIFluid('phases', 'WO', ...
                           'mu',    [1, 5]*centi*poise, ...
                           'rho',   [1000, 700]*kilogram/meter^3, ...
                           'n',     [2, 2]);

% Once setup, we can modify the fluid structure and replace some of the
% functions. Here we implement a function to give constant oil
% compressibilty.
c        = 0.001/barsa;
p_ref    = 300*barsa;
fluid.bO = @(p, varargin) exp((p - p_ref)*c);
disp(fluid)


%% Construct reservoir simulator class
% We create the simulation model using the TwoPhaseOilWaterModel. It is
% important to make sure gravity is turned on before creating the model
% otherwise gravity may not be included.
gravity reset on
model = TwoPhaseOilWaterModel(G, rock, fluid);

% Once the model is created we can inspect the fluid properties
% graphically.
inspectFluidModel(model,'field','Densities')

%% Define initial state
% We setup the initial state (initial saturation and pressure in each cell)
% by specifying the depth of the oil water contact aand the pressure and
% the datup depth. The state will be setup with water at the bottom and oil
% on the top.
depthOW = 85*meter;
depthD  = 10*meter;
region = getInitializationRegionsBlackOil(model, depthOW, ...
            'datum_depth', depthD, 'datum_pressure', p_ref);
state0 = initStateBlackOilAD(model, region);

% Visualise
clf;
plotCellData(G, state0.s(:,1),'EdgeAlpha',.1), colormap(flipud(winter));
plotWell(G,W,'color','k')
patch([-50 1050 1050 -50],[-50 -50 1050 1050],depthOW*ones(1,4), ...
    ones(1,4), 'FaceColor',[.6 .6 1],'EdgeColor','r','LineWidth',1);
view(-50, 50), axis tight off


%% Simulation schedule and set solver parameters
% Compute the timesteps. rampupTimesteps creates timesteps that ramp up 
% geometrically up to a specificed time.
nstep   = 25;
timesteps = rampupTimesteps(simTime,simTime/(nstep + 1),5);

% Set up the schedule containing both the wells and the timestep
schedule = simpleSchedule(timesteps, 'W', W);
 
% The default tolerances can be seen by inspecting the model object. For
% our particular setup there are convergence issues when using the
% defaults so we need to tighten some of the tolerances.
model.drsMaxRel = inf;
model.dpMaxRel  = .1;
model.dsMaxAbs  = .1;

% We can automatically determine the best solver to use for the model. Here
% the linear solver is chosen to be the AMGCL_CPRSolverAD. This is a
% compiled solver which will run the simulation faster than if we were
% using a standard backslash method. If AMGCL is not installed and compiled
% on your system, mrst will try to do this automatically. 
solver = getNonLinearSolver(model);

%% Run and monitor simulation
% We setup a function which does some plotting after every timestep. This
% is passed in the simulateScheduleAD function as the 'afterStepFn', which 
% is a function which runs after every step.

wsargs = {'field', 'wcut', 'SelectedWells', [1,2,3]};
fn = getPlotAfterStep(state0, model, schedule,'view',[50 50], ...
                     'field','s:1','wells',W,'plotWellSolArgs',wsargs);

% simulateScheduleAD will run through the timesteps setup in the schedule
% and produce results for each timestep. wellSols contains results for the
% wells. States contains the primary variables at each timesteps and report
% contains output from the solver. 
[wellSols, states, report] = ...
   simulateScheduleAD(state0, model, schedule, ...
                       'NonLinearSolver', solver, 'afterStepFn',fn);

%% Copyright notice

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

