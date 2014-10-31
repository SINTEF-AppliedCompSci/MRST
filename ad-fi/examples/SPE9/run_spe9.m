%% Ninth Comparative Solution Project
% This example runs the model from Killough, J. E. 1995. Ninth SPE
% comparative solution project: A reexamination of black-oil simulation. In
% SPE Reservoir Simulation Symposium,  12-15 February 1995, San Antonio,
% Texas. SPE 29110-MS, doi: 10.2118/29110-MS
try
   require ad-fi deckformat
catch %#ok<CTCH>
   mrstModule add ad-fi deckformat
end

mrstVerbose true

%% Read and process file.
current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'SPE9.DATA');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

fluid = initDeckADIFluid(deck);
gravity on

%% Initial solution
p0  = deck.SOLUTION.PRESSURE;
sw0 = deck.SOLUTION.SWAT;
sg0 = deck.SOLUTION.SGAS;
s0  = [sw0, 1-sw0-sg0, sg0];
rs0 = deck.SOLUTION.RS;
rv0 = 0;

state = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);   clear k p0 s0 rs0

schedule = deck.SCHEDULE;
system = initADISystem(deck, G, rock, fluid, 'cpr', true);
% switch off individual well solves:
system.stepOptions.solveWellEqs = false;
% use new cpr based on dnr
system.nonlinear.cprBlockInvert = false;
% convergence is overall better for quite strict limits on update
system.stepOptions.drsMax = .2;
system.stepOptions.dpMax  = .2;
system.stepOptions.dsMax  = .05;
% gmres tol needs to be quite strict
system.nonlinear.cprRelTol = 1e-3;
system.pscale = 1/(200*barsa);

%% Plot reservoir and wells
figure, 
plotCellData(G,log10(convertTo(rock.perm(:,1),milli*darcy))); view(3); 
h = colorbar('horiz');
set(h,'XTickLabel', num2str(10.^(get(h,'XTick')')), ...
   'YTick', .5, 'YTickLabel','[mD]','Position',[.13 .05 .77 .03]);
 W = processWellsLocal(G,rock, schedule.control(1));
plotWell(G,W(1),'color','b','linewidth',3,'fontsize',12);
plotWell(G,W(2:end),'color','k','linewidth',3,'fontsize',12);
axis tight, view(40,55);

%% Run schedule
timer = tic;
[wellSols, states, iter] = runScheduleADI(state, G, rock, system, schedule);
toc(timer)

%% Plot solutions
T     = convertTo(cumsum(schedule.step.val), day);
[qWs, qOs, qGs, bhp] = wellSolToVector(wellSols);
injInx = 1;
prdInx = 2:26;
figure(1)
plot(T, qWs(:, injInx)*day, '-*'), xlabel('Years'),title('Water injection rate')
figure(2)
plot(T, bhp(:, prdInx)/barsa, '-*'), xlabel('Years'),title('BHPs')
figure(3)
plot(T, qOs(:, prdInx)*day, '-*'), xlabel('Years'),title('Oil production rates')
figure(4)
plot(T, qGs(:, prdInx)*day, '-*'), xlabel('Years'),title('Gas production rates')







