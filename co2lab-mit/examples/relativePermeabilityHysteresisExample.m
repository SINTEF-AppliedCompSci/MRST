%% Black Oil Test: relative permeability hysteresis
%
% 3D mesh with 2 fluids (water and CO2)
%
% Objective:
% Run a simulation with 1 region, no dissolution, with and without
% hysteresis. Assuming that the formation is strongly water-wet, kr
% hysteresis is only considered for the CO2 phase. Accounting for
% hysteresis is expected to greatly reduce the volume of free CO2 at the
% top of the injection layer, since residual trappping occurs as the
% plume migrates upwards. PUNQ-S3 model, see Juanes et al., WRR 2006.

%% Initial settings
% Modules, etc
clear, clc, close all
mrstModule add ad-props deckformat mrst-gui ad-core ad-blackoil linearsolvers
gravity reset on;
mrstVerbose off

% Model options
krhyst = 0;         % use relative permeability hysteresis
use_mex = 0;
use_cpr = 0;
fine_tstep = 0;     % use finer timesteps than deck

fn = fullfile(getDatasetPath('co2labmit'), 'punq-s3', 'CASE2.DATA');

if krhyst == 1
    folderName = 'case2_hyst';
else
    folderName = 'case1_nohyst';
end

% Plotting options
plot_grid  = 0;
plot_permX = 0;
plot_kr    = 0;
plot_pvm   = 0;
plot_p0    = 0;


%% ------------- Input deck and generate grid, rock and fluid -------------
% We use the same input deck (CASE2.DATA) and properties that was used in 
% Juanes et al., WRR (2006), since MRST is able to import ECLIPSE decks.
deck = convertDeckUnits(readEclipseDeck(fn));
G = initEclipseGrid(deck);
shiftZup = 840 - min(G.nodes.coords(:, 3));
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + shiftZup;
G = computeGeometry(G);
if plot_grid == 1
    plotGrid(G); set(gca,'YDir','reverse');  view([-40 66])
end

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
%rock.poro = rock.poro + (0.2 - mean(rock.poro)); % mean poro = 0.2
if plot_permX == 1
    plotCellData(G, rock.perm(:,1)/(milli*darcy)); set(gca,'YDir','reverse');
    view([-40 66]); colormap(jet); c=colorbar; caxis([0.5 1000])
end

fluid = initDeckADIFluid(deck);
%fluid.bO = [];                                 % needed for mex
rock.regions.saturation = ones(G.cells.num, 1); % 1 saturation region
if krhyst == 1
    fluid.krHyst = 2;   % imbibition regions where hysteresis is active  
else
    % no hysteresis
    fluid = rmfield(fluid, 'krHyst');
end

% --------------------- Plot hysteretic curve ----------------------------
if plot_kr == 1
    sgd = linspace(fluid.krPts.g(1,2),fluid.krPts.g(1,3),50)';
    krd = fluid.krG{1}(sgd);
    sgi = linspace(fluid.krPts.g(2,2),fluid.krPts.g(2,3),50)'; 
    kri = fluid.krG{2}(sgi);
    sgs = linspace(0.28,0.4,50)';
    krs = fluid.krGi{1}(sgs,repelem(0.4,50,1));
    figure(2)
    plot(sgd,krd,'-r','linewidth', 1.5, ...
         'DisplayName', '$k_\mathrm{r,g}^\mathrm{d}$'), hold on
    plot(sgi,kri,'--r','linewidth', 1.5, ...
         'DisplayName', '$k_\mathrm{r,g}^\mathrm{i}$')
    plot(sgs,krs,'-.r','linewidth', 1, ...
         'DisplayName', '$k_\mathrm{r,g}^\mathrm{s}$'), hold off
    grid on
    xlabel('$S_\mathrm{g}$ [-]', 'interpreter', 'latex')
    ylabel('$k_\mathrm{r,g}$ [-]', 'interpreter', 'latex')
    title('Hysteretic gas relative permeability (Berea sst., gas-water, from Oak, 1990)')
    legend('interpreter', 'latex', 'fontsize', 12)
    xlim([0 1]); ylim([0 1])
end


%% ----------------------- Model and ancceleration ------------------------
model = selectModelFromDeck(G, rock, fluid, deck); 
if use_mex == 1
    model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);
end
model = model.validateModel();  % Set up the state function groups
if use_cpr == 1
    nls = getNonLinearSolver(model, 'TimestepStrategy', 'iteration', ...
                             'useCPR', true);
    nls.LinearSolver = AMGCL_CPRSolverAD();
else
    nls = getNonLinearSolver(model, 'TimestepStrategy', 'iteration');
end
nls.useLinesearch = true;
nls.maxIterations = 15; 
nls.maxTimestepCuts = 12;
nls.acceptanceFactor = 2;

% Pore volume multipliers
% Note that here, we multiply the pore volume of the boundary cells by
% 1000, following Juanes et al. (2006)
idAct = find(deck.GRID.ACTNUM);
idG = (1:prod(deck.GRID.cartDims))';
idG_mult = idG(deck.GRID.PORV==1000);
cellsb = ismember(idAct, idG_mult);
if plot_pvm == 1
    plotGrid(G, 'facecolor', 'none'); 
    plotGrid(G,cellsb);
    set(gca,'YDir','reverse');  view([-40 66])
end
model.operators.pv(cellsb) = model.operators.pv(cellsb)*1e3;
model.FlowPropertyFunctions = ...
model.FlowPropertyFunctions.setStateFunction('RelativePermeability', ...
                                             HystereticRelativePermeability(model));
% Specify Outputs (these can also be computed after the simulation run via model.getProp)
%model.OutputStateFunctions = {'ComponentTotalMass','RelativePermeability', 'Density', ...
%                              'Mobility', 'PhasePressures', 'PoreVolume', 'ShrinkageFactors', 'Viscosity'};


%% -------------------------- Initialization -------------------------------
g = norm(gravity);
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
p_r = g*fluid.rhoWS*z_0;   
equil  = ode23(@(z,p) g .* fluid.bW(p)*fluid.rhoWS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
s0  = repmat([1, 0], [G.cells.num, 1]);  % fully saturated in water --> [1 0]
rs0 = 0;                                 % immiscible 
rv0 = 0;                                 % dry gas
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

if plot_p0 == 1
    plotCellData(G, p0/barsa); 
    set(gca,'YDir','reverse');  view([-40 66])
    colormap(turbo), c = colorbar; c.Label.String = '$p_0$ [bar]'; 
    c.Label.Interpreter = 'latex'; c.Label.FontSize = 12;
end


%% ------- Deck schedule into a MRST schedule by parsing the wells --------
% Juanes et al., WRR (2006) inject at a rate of 18 rm3/day/well over 10y, and run the 
% simulation for 500y (Fig. 5)

% Read schedule
schedule = convertDeckScheduleToMRST(model, deck);

% needed schedule updates
% Update well rate to corresponding value at surface conditions, since running
% the case with 'resv' type well will give an error.
resv = 18/day; % this is the value in schedule (m3/s per well)
pavg = mean(p0([schedule.control(1).W.cells]));
rate = resv*fluid.bG(pavg); % rm3 to sm3
[schedule.control(1).W.val] = deal(rate);
[schedule.control(2).W.val] = deal(0);          % control 2 no injection
[schedule.control(1).W.type] = deal('rate');    % supported types are 'rate' or 'bhp'
[schedule.control(2).W.type] = deal('rate');     
nwell = numel(schedule.control(1).W);
for n=1:nwell                                   % add the p limit as in the paper (not reached)
    schedule.control(1).W(n).lims.bhp = 160*barsa;
end

% Modify timesteps?
if fine_tstep == 1
    repTimes = [0.05:0.05:6 6.1:.1:10 10.5:0.5:20 21:50 52:2:100 ...
                105:5:200 210:10:500]'*year;
    timesteps = [repTimes(1); diff(repTimes)];
    schedule.step.val = timesteps;
    % schedule.step.control = ones(numel(timesteps), 1)+1;
    % schedule.step.control(1:20) = 1;
    schedule.step.control = ones(numel(timesteps), 1);
    schedule.step.control(161:end) = 2;
end

% Plot permeability and wells (compare to Fig. 2 in Juanes et al., WRR
% 2006). Note changes in the colormap!
if plot_permX == 1
    figure;
    plotCellData(G, convertTo(rock.perm(:,1), milli*darcy), 'FaceAlpha', 0.95, 'EdgeAlpha', 0.3, 'EdgeColor', 'k');
    plotWell(G, schedule.control(1).W, 'radius',.5); % Pick the only well control present
    set(gca,'YDir','reverse')
    title('Permeability (mD)')
    axis tight, view(-40, 70), c = colorbar('SouthOutside'); caxis([0.5 1000]);
    colormap(jet)
end


%% ------------------------------ Simulation -------------------------------
N=1;
maxNumCompThreads(N)
if use_cpr == 1
    nls.LinearSolver.amgcl_setup.nthreads = N; 
end
problem = packSimulationProblem(state0, model, schedule, folderName,...
        'Name', folderName, 'NonLinearSolver', nls);
[ok, status] = simulatePackedProblem(problem);
[wellSols, states, report] = getPackedSimulatorOutput(problem);  


%% ------------------------------- Results --------------------------------
% approximate ECLIPSE colormap 
colors = [
    0   0   1       % Blue
    0   1   1       % Cyan
    0   1   0       % Green
    1   1   0       % Yellow
    1   0   0       % Red
    ];
cmap1 = interp1([0 1], colors(1:2, :), linspace(0,1,30)); 
cmap2 = interp1([0 1], colors(2:3, :), linspace(0,1,25)); 
cmap4 = interp1([0 1], colors(3:4, :), linspace(0,1,25)); 
cmap5 = interp1([0 1], colors(4:5, :), linspace(0,1,20)); 
cmap = [cmap1; cmap2; cmap4; cmap5];

% Comparison figure
figure(11)
plotCellData(G, states{end}.s(:, 2), ...
             'FaceAlpha', 1, 'EdgeAlpha', 0.3, 'EdgeColor', 'k');
plotWell(G, schedule.control(1).W, 'radius',.5); % Pick the only well control present
set(gca,'YDir','reverse')
%title(['CO_2 saturation, t=' num2str(round(sum(schedule.step.val)/year)) 'y']); 
colormap(cmap)
axis tight off, view(-40, 66), c = colorbar('SouthOutside'); caxis([0 1]);
c.Ticks = [0 0.3 0.5 0.8 1];

% Other results
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;  
    %pc = model.getProp(states{n}, 'CapillaryPressure');
    %states{n}.pc = pc{1,2};
end
figure(10)
plotToolbar(G, states); view(-40, 70); set(gca,'YDir','reverse')

% kr Hyst plots 
% ns = 50;
% s = linspace(0, 0.7, ns);
% kr = zeros(ns, ns);
% kr_d = zeros(ns, ns);
% for i = 1:ns
%     for j = 1:ns
%         sg = s(i);
%         sg = initVariablesADI(sg);
%         sgmax = s(j);
%         sgmax = max(sgmax, sg);
%        
%         kg = fluid.krGi{end}(sg, sgmax);
%         kr(i, j) = value(kg);
%         kr_d(i, j) = kg.jac{1}(1, 1);
%     end
% end
% figure;
% plot(s, kr, '.-');
% title('kr_g')
% 
% figure;
% stairs(s, kr_d);
% title('d kr_g / d sg')