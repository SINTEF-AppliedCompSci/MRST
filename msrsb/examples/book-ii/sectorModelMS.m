%% MsRSB as a CPR preconditioner
% This example uses the simulation model from the workflow example
% discussed in Section 15.6.4 of the MRST book (upscalingExample3.m in the
% book module) to demonstrates how we can use MsRSB in MRST's fully
% implicit framework. That is, we use the iterative method as the elliptic
% solver in a CPR preconditioner.

mrstModule add ad-core ad-blackoil ad-props mrst-gui msrsb coarsegrid linearsolvers
% Set a single thread for all solvers to leave parallelism out of the
% discussion. You can increase this number to reduce the runtime.
t0 = maxNumCompThreads(1);
%% Build simulation model
% The simulation model is built manually. Here, we only give the necessary
% statements and refer to the original example for discussion and plots.
rng(0);
[xmax,ymax, n]  = deal(1000*meter, 1000*meter, 30);
[x, y]  = meshgrid(linspace(0,xmax,n+1), linspace(0,ymax,n+1));
[x, y]  = deal(x',y');
dome    = 1-exp(sqrt((x - xmax/2).^2 + (y - ymax/2).^2)*1e-3);
[xn,yn] = deal(pi*x/xmax,pi*y/ymax);
perturb = sin(5*xn) + .5*sin(4*xn+6*yn) + cos(.25*xn*yn./pi^2) + cos(3*yn);
perturb = perturb/3.5;
[h, hr] = deal(8,1);
zt      = 50 + h*perturb + rand(size(x))*hr - 20*dome;
zb      = zt + 30;
zmb     = min(zb + 4 + 0.01*x - 0.020*y + hr*rand(size(x)), zb);
zmt     = max(zb -15 + 0.01*x - 0.025*y + hr*rand(size(x)), zt);

horizons = {struct('x', x, 'y', y, 'z', zt), ...
            struct('x', x, 'y', y, 'z', zmt), ...
            struct('x', x, 'y', y, 'z', zmb), ...
            struct('x', x, 'y', y, 'z', zb)};
grdecl   = convertHorizonsToGrid(horizons, 'dims', [40 40], 'layers', [3 6 3]);
[X,Y,Z]  = buildCornerPtNodes(grdecl);
i=47:80; Z(i,:,:) = Z(i,:,:) + .022*min(0,Y(i,:,:)-550);
j= 1:30; Z(:,j,:) = Z(:,j,:) + .021*min(0,X(:,j,:)-400);
j=57:80; Z(:,j,:) = Z(:,j,:) + .023*min(0,X(:,j,:)-750);
grdecl.ZCORN = Z(:);

G = computeGeometry(processGRDECL(grdecl));

% Petrophysics
rng(357371);
[K,L] = logNormLayers(G.cartDims, [100 400 10 50]*milli*darcy);
K     = K(G.cells.indexMap);
perm  = [K, K, 0.1*K];
rock  = makeRock(G, perm, 0.3);

% Define wells
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
W = verticalWell(W, G, rock, G.cartDims(1)-5, offset, [],...
                'Name', 'I1', 'comp_i', [1 0], ...
                'Val', injRate, 'Type', 'rate', 'refDepth', 50);

% Three-phase template model with constant oil compressibility
fluid = initSimpleADIFluid('mu',    [1, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);
c        = 0.001/barsa;
p_ref    = 300*barsa;
fluid.bO = @(p, varargin) exp((p - p_ref)*c);

% Construct reservoir model and initial state
gravity reset on
model  = TwoPhaseOilWaterModel(G, rock, fluid);
region = getInitializationRegionsBlackOil(model, 85*meter, ...
            'datum_depth', 10*meter, 'datum_pressure', p_ref);
state0 = initStateBlackOilAD(model, region);

% Define simulation schedule and set solver parameters
nstep      = 25;
startSteps = repmat((simTime/(nstep + 1))/5, 5, 1);
restSteps  = repmat(simTime/(nstep + 1), nstep, 1);
timesteps  = [startSteps; restSteps];
schedule   = simpleSchedule(timesteps, 'W', W);

% Tighten tolerences
model.drsMaxRel = inf;
model.dpMaxRel  = .1;
model.dsMaxAbs  = .1;

% Set up CPR preconditioner
pressureSolver = AMGCLSolverAD('tolerance', 1e-4, ...
                               'coarsening', 'aggregation', ...
                               'solver', 'gmres', ...
                               'reuseMode', 2); % Reuse setup inside CPR
try
    % Check if solver is working
    pressureSolver.solveLinearSystem(speye(3), rand(3, 1));
catch
    % Fall back to Matlab's builtin direct solver
    pressureSolver = BackslashSolverAD();
end
linsolve = CPRSolverAD('ellipticSolver', pressureSolver, 'tolerance', 1e-3);

%%
figure(1)
K = convertTo(rock.perm(:,1),milli*darcy);
plotCellData(G, log10(K),'EdgeAlpha',.1);
mrstColorbar(K,'east',true);
plotWell(G,W,'FontSize',12);
view(-75,50); axis tight off

%% Simulate base case
[wellSols, states, report] = ...
   simulateScheduleAD(state0, model, schedule, 'LinearSolver', linsolve);

%% Simulate with multiscale solver
% Make coarse grid
p  = partitionUI(G, [8 8 3]);
p  = processPartition(G, compressPartition(p));
CG = coarsenGeometry(generateCoarseGrid(G, p));
CG = storeInteractionRegion(CG);
CG = setupMexInteractionMapping(CG);
figure(1); 
plotFaces(CG, 1:CG.faces.num,'FaceColor','none','EdgeColor','k');

% Set up multiscale solver as CPR preconditioner
msSolver = MultiscaleVolumeSolverAD(CG, 'tolerance', 1e-4, ...
            'maxIterations', 1, 'useGMRES', false, 'verbose', false, ...
            'getSmoother', getSmootherFunction('type', 'ilu0', 'iterations', 1));

linsolve = CPRSolverAD('ellipticSolver', msSolver, 'tolerance', 1e-3);
% Solve problem
[wellSolsMS, statesMS, reportMS] = ...
   simulateScheduleAD(state0, model, schedule, 'LinearSolver', linsolve);

%% Plot comparison
plotWellSols({wellSols, wellSolsMS}, cumsum(schedule.step.val), ...
    'datasetnames',{'amg','msrsb'}, 'field','qOs');

%% Plot time consumption
% Get total runtime
timeMS   = cellfun(@(x) x.WallTime, reportMS.ControlstepReports);
timeBase = cellfun(@(x) x.WallTime, report.ControlstepReports);

% Get detailed timining
solveMS   = getReportTimings(reportMS);
solveBase = getReportTimings(report);

solveTiming = [vertcat(solveMS.LinearSolve), vertcat(solveBase.LinearSolve)];
solveTotal = [timeMS, timeBase];

figure, hold all
bar(vertcat(solveBase.LinearSolve))
bar(vertcat(solveMS.LinearSolve),'BarWidth',.3)
legend('amg','msrsb'); set(gca,'FontSize',12);
xlabel('Time step'), ylabel('Linear solver time [s]');

figure, hold all
bar(timeBase)
bar(timeMS,'BarWidth',.3)
legend('amg','msrsb'); set(gca,'FontSize',12);
xlabel('Time step'), ylabel('Runtime [s]');

%% Solve using sequential multiscale solver
% We can also use the multiscale method in a sequential framework. We build
% a sequential model and set parameters to get exact mass conservation for
% a single pass of pressure and transport with the same convergence
% criterion as the fully-implicit solver
mrstModule add sequential
seqmodel = getSequentialModelFromFI(model);

seqmodel.pressureModel.incTolPressure = 1e-5;
seqmodel.transportModel.conserveOil = true;
seqmodel.transportModel.conserveWater = true;
seqmodel.transportModel.useCNVConvergence = true;
seqmodel = addMultiscaleSolverComp(seqmodel, CG, 'maxIterations', 50,...
                                               'useGMRES', true, ...
                                               'tolerance', 0.01);

% Run sequentially implicit
wsSIms = simulateScheduleAD(state0, seqmodel, schedule);

%% Run sequentially fully implicit
% Here, we configure the algorithm to only check the volume discrepancy and
% not the fully implicit residual or the increments in saturations from one
% transport step to the next within the outer iteration. The first line
% tells the simulator that the time step should be iterated upon.
seqmodel.stepFunctionIsLinear        = false;
seqmodel.outerCheckParentConvergence = false;
seqmodel.volumeDiscrepancyTolerance  = 3e-3;
seqmodel.incTolSaturation            = inf;
wsSFIms = simulateScheduleAD(state0, seqmodel, schedule);

%% Compare the solutions
plotWellSols({wellSols,wsSFIms,wsSIms}, cumsum(schedule.step.val), ...
    'datasetnames',{'FI','SFI','SI'},'field','qOs');
%% Reset the threads back to default
maxNumCompThreads(t0);
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
