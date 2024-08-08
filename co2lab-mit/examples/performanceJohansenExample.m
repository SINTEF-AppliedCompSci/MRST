%% Conceptual Example: CO2 Injection in the Johansen formation, North Sea
% Note: This conceptual example builds on the showJohansen.m example from 
%       the MRST book. The fluid properties were selected for the purposes
%       of illustrating MRST usage in a realistic case only, and are not
%       meant to be representative of the actual properties of the Johansen 
%       formation.

% Setup
mrstModule add ad-props ad-blackoil deckformat ad-core mrst-gui linearsolvers co2lab-mit
mrstVerbose on

% Data input
directory = fullfile(getDatasetPath('co2labmit'), 'example_Johansen/');
fn = 'example_Johansen.DATA';

% Options
krhyst = true;              % set to true for kr hysteresis
useAcceleration = true;     % use compiled linear solver and AD backend
injTime = 30*year;          % s
simTime = 500*year;         % s
mrate = 1*mega*tonne/year;  % kg/s
plot_figs = false;

% Output
if krhyst == true
    folderName = 'example_Johansen_hyst';
else
    folderName = 'example_Johansen_nohyst';
end

if useAcceleration
    caseName = 'accelerated';
else
    caseName = 'baseline';
end

%% Grid
% We start by reading the model from a file in the Eclipse format (GRDECL),
% picking the sector model with five vertical layers in the Johansen
% formation and with five shale layers above and one below.
dpath = getDatasetPath('johansen');
sector = fullfile(dpath, 'NPD5');
filename = [sector, '.grdecl'];
grdecl   = readGRDECL(filename);
G = computeGeometry(processGRDECL(grdecl));

if plot_figs == true
    % Plot grid (all)
    grdecl_all = grdecl;
    grdecl_all.ACTNUM = ones(numel(grdecl.ACTNUM), 1);
    G_all = processGRDECL(grdecl_all);

    clf, 
    h = plotGrid(G_all,'FaceColor',.9*[1 1 1], 'FaceAlpha',.5,'EdgeColor',[0.2 0.2 0.2], 'EdgeAlpha', 0.5);
    plotFaces(G_all,find(G_all.faces.tag>0),'FaceColor','r','EdgeColor','none');
    axis tight off; view(-145,60);
    set(gca,'ZDir','normal'), camlight headlight, set(gca,'Zdir','reverse');
    
    % Plot depth map
    clf,
    plotCellData(G,G.cells.centroids(:,3),'EdgeColor','k','EdgeAlpha',0.1);
    view(3), axis tight off, view(-145,60), zoom(1.2)
    set(gca,'Clipping','off'), colormap(turbo)
    c = colorbar('southoutside'); c.Label.String = 'Depth [m]'; 
    c.Label.Interpreter = 'latex'; c.FontSize = 14;
    c.Label.FontSize = 18;
     
end


%% Porosity and Permeability
% The porosity data are given with one value for each cell in the model. We
% read all values and then pick only the values corresponding to active
% cells in the model.
% The permeability is given as a scalar field (Kx) similarly as the
% porosity. The tensor is given as K = diag(Kx, Kx, 0.1Kx).

% Assign to rock object
p = reshape(load([sector, '_Porosity.txt'])', prod(G.cartDims), []);
poro = p(G.cells.indexMap); clear p
K = reshape(load([sector, '_Permeability.txt']')', prod(G.cartDims), []);
perm = bsxfun(@times, [1 1 0.1], K(G.cells.indexMap)).*milli*darcy; clear K;
rock = makeRock(G, perm, poro);

if plot_figs == true
    % Visualize porosity
    clf
    subplot(1,2,1)
    p = reshape(load([sector, '_Porosity.txt'])', prod(G.cartDims), []);
    plotCellData(G, rock.poro,'EdgeColor','k','EdgeAlpha',0.1);
    colorbarHist(rock.poro,[0.09 0.31],'West',50); view(-45,15), axis tight off, zoom(1.2)
    set(gca,'Clipping','off'); colormap(turbo)
    subplot(1,2,2)
    view(-15,40)
    plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
    plotCellData(G, poro, poro>0.1, 'EdgeColor','k','EdgeAlpha',0.1);
    colorbarHist(poro(poro>.1),[0.09 0.31],'West',50);
    axis tight off, zoom(1.45), camdolly(0,.2,0);
    set(gca,'Clipping','off')
    
    % Visualize permeability
    clf
    subplot(1,2,1)
    p = log10(rock.perm(:,1)/(milli*darcy));
    plotCellData(G,p,'EdgeColor','k','EdgeAlpha',0.1);
    view(-45,15), axis off,
    set(gca,'Clipping','off'); colormap(turbo)
    cs = [-2.5 -2 -1 0 1 2 3];
    clim([min(cs), max(cs)]);
    h = colorbarHist(p,[min(cs), max(cs)],'South',25);
    h.Ticks = cs; h.Label.String = 'log$_{10} (k$ [mD])'; 
    h.Label.Interpreter = 'latex'; h.FontSize = 14;
    h.Label.FontSize = 18;
    subplot(1,2,2)
    idx = p>1;
    plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
    plotCellData(G, p, idx, 'EdgeColor','k', 'EdgeAlpha', 0.1);
    view(-20,35), axis off, 
    set(gca,'Clipping','off')
    cs = [1 2 3];
    clim([min(cs), max(cs)]);
    h = colorbarHist(p(idx), [min(cs), max(cs)], 'South', 25);
    h.Ticks = cs; h.Label.String = 'log$_{10} (k$ [mD])'; 
    h.Label.Interpreter = 'latex'; h.FontSize = 14; 
    h.Label.FontSize = 18;
end


%% Fluid && Regions
deck = convertDeckUnits(readEclipseDeck([directory, fn]));
fluid = initDeckADIFluid(deck);

nreg = 2;
id_shale = log10(rock.perm(:,1)/(milli*darcy)) < -0.5;
rock.regions.saturation = ones(G.cells.num, 1);
rock.regions.saturation(id_shale) = 2;
rock.regions.imbibition = rock.regions.saturation + nreg;

% Hysteresis
if krhyst == 1
    fluid.krHyst = 3;   % imbibition regions where hysteresis is active  
else
    % no hysteresis
    fluid = rmfield(fluid, 'krHyst');
end

if plot_figs == true
    % Visualize regions
    id_reg = rock.regions.saturation;
    clf
    ax(1) = subplot(1,2,1);
    plotCellData(G,id_reg,'EdgeColor','k','EdgeAlpha',0.1);
    view(-45,15), axis tight off, zoom(1.2)
    set(gca,'Clipping','off'); 
    cs = [1 2];
    caxis([min(cs), max(cs)]);
    h = colorbar; h.Ticks = [1 2]; h.Label.String = 'Region number';
    ax(2) = subplot(1,2,2);
    idx = rock.regions.saturation == 1;
    plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
    plotCellData(G, id_reg, idx, 'EdgeColor','k', 'EdgeAlpha', 0.1);
    view(-20,35), axis tight off, zoom(1.45), camdolly(0,.12,0);
    set(gca,'Clipping','off'); 
    colormap(ax(1), turbo)
    colormap(ax(2), flipud(turbo(2)))
end


%% Initialization
clear p
gravity reset on
g = norm(gravity);
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
p_r = g*fluid.rhoOS*z_0;   
equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS, [z_0, z_max], p_r);
p0  = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
s0  = repmat([1, 0], [G.cells.num, 1]);  % fully saturated in water --> [1 0]
rs0 = zeros(G.cells.num, 1);             % immiscible 
rv0 = 0;                                 % dry gas
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

% Plot
if plot_figs == true
    clf
    plotCellData(G, p0/barsa, 'EdgeColor','k','EdgeAlpha',0.1); 
    view(-45,15), axis tight off
    colormap(turbo)
    c = colorbar; c.Label.String = '$p_0$ [bar]'; 
    c.Label.Interpreter = 'latex'; c.Label.FontSize = 12;
end


%% Wells
injrate = mrate/fluid.rhoGS;  % m3/s (surface volume)

% Create well
w = load([sector, '_Well.txt']);
W = verticalWell([], G, rock,  w(1,1), w(1,2), w(1,3):w(1,4),  ...
                 'InnerProduct', 'ip_tpf', 'Radius', 0.1, 'name', 'I');
% Specify rate
W.type = 'rate';
W.val = injrate;
W.compi = [0 1];
W.refDepth = G.cells.centroids(W.cells(1), G.griddim);

% add to previous plot
if plot_figs == true
    plotWell(G,W,'height',500,'color','r');
end


%% Model
model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, 'water', false);

% Acceleration
nls = getNonLinearSolver(model, 'TimestepStrategy', 'iteration');
if useAcceleration
    model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, 'rowMajor', true);
else
    nls.LinearSolver = BackslashSolverAD(); % Override automatic choice
end
model = model.validateModel();  % Set up the state function groups
nls.maxIterations = 12; 

% BCs and schedule
bc = [];

reportTimes = [1*hour, [1, 2, 7, 14, 30, 60, 120, 180, 365.25]*day, ...
              [2:.5:5 6:1:50 52:2:100 105:5:500]*year];
timesteps = [reportTimes(1) diff(reportTimes)];
assert(sum(timesteps)==simTime, 'sum of timesteps must equal simTime')

schedule_inj = simpleSchedule(timesteps, 'W', W, 'bc', bc);                 % Simple schedule, this would inject for the total simTime
tmp = cell(2, 1);                                                           % create 2 schedules
schedule = struct('step', schedule_inj.step);                               % timesteps and wells for each timestep
schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);                 % add 2 fields for 2 wells
schedule.control(1).W = W;                                                  % field 1 used during injection
schedule.control(2).W = W;                                                  % nr of wells must be the same for the entire simulation
schedule.control(2).W.val = 0;                                              % field 2 rate 0 (after injection)
schedule.control(1).bc = bc;
schedule.control(2).bc = bc;
schedule.step.control(cumsum(schedule.step.val) > injTime) = 2;             % set timesteps after injTime to use Well field 2


%% Run
close all
problem = packSimulationProblem(state0, model, schedule, folderName, 'Name', caseName, 'NonLinearSolver', nls);
[ok, status] = simulatePackedProblem(problem);
[wellSols, states, report] = getPackedSimulatorOutput(problem); 


%% Results
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;  
end

% Overview
p = log10(rock.perm(:,1)/(milli*darcy));
idx = p > 1;
clf
plotToolbar(G, states, idx); view(-40, 70);
view(-45,15), axis tight off
colormap(turbo)

% Summary figures
figure(2)
plotCellData(G, states{42}.dp/barsa, idx, 'EdgeColor', 'k', 'EdgeAlpha', 0.1);
plotWell(G,W,'height',100,'color','k');
view(3), axis tight off, view(-145,60), zoom(1.2)
set(gca,'Clipping','off'), colormap(turbo)
c = colorbar('southoutside'); c.Label.String = '$\Delta p$ [bar]'; 
c.Label.Interpreter = 'latex'; c.FontSize = 14;
c.Label.FontSize = 18;

figure(3)
plotCellData(G, states{end}.s(:,2), idx, 'EdgeColor', 'k', 'EdgeAlpha', 0.1);
view(3), axis tight off, view(-145,60), zoom(1.2)
set(gca,'Clipping','off'), colormap(turbo)
caxis([0 0.8])
c = colorbar('southoutside'); c.Label.String = '$S_g$ [-]'; 
c.Label.Interpreter = 'latex'; c.FontSize = 14;
c.Label.FontSize = 18;

%% Performance
% We compare the performance of the base case (direct solver, standard
% sparse backend) and the accelerated version (iterative solver,
% mex-accelerated backend).
problem_accelerated = packSimulationProblem(state0, model, schedule, folderName, 'Name', 'accelerated');
problem_base = packSimulationProblem(state0, model, schedule, folderName, 'Name', 'baseline');

[~, ~, reports_accelerated] = getPackedSimulatorOutput(problem_accelerated);
[~, ~, reports_base] = getPackedSimulatorOutput(problem_base);
%%
allreports = {reports_base, reports_accelerated};
allnames = {'Direct solver, standard backend', 'AMGCL solver, diagonal AD backend'};
timings = cell(2, 1);
active = true(2, 1);
data = zeros(2, 3);
for i = 1:2
    r = allreports{i};
    if isempty(r)
        active(i) = false;
        continue
    end
    name = allnames{i};
    t = getReportTimings(r, 'total', true);
    lsolve = t.LinearSolve + t.LinearSolvePrep;
    fprintf('\n%s:\n', name)
    fprintf('----------------------------------\n')
    fprintf('Linear solver: %5.4f seconds (%2.4f s / Newton)\n', lsolve, lsolve./t.Iterations)
    fprintf('Eqn. Assembly: %5.4f seconds (%2.4f s / Assembly)\n', t.Assembly, t.Assembly/t.NumberOfAssemblies)
    fprintf('Total runtime: %5.4f seconds (%2.4f s / Newton) \n', t.Total, t.Total/t.Iterations)
    data(i, 1) = lsolve;
    data(i, 2) = t.Assembly;
    data(i, 3) = t.Total;
    timings{i} = t;
end
