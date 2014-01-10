
mrstModule add ad-fi deckformat spe10 mex mrst-gui

gravity on
% Layer indices
layers = 1:5;

% Capping of pore volume to avoid messing with actnum
minporo = 0.01;

% Read and process file.
current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'SPE10_OW.DATA');

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);


[G, W, rock] = SPE10_setup(layers);
rock.poro(rock.poro < minporo) = minporo;
fluid = initDeckADIFluid(deck);

state0 = initResSol(G, 6000*psia, [.01 .99 0]);
state0.rs = 0*state0.pressure;

for i = 1:4;
    W(i).val = 4000*psia;
    W(i).lims.bhp = -inf;
    W(i).lims.rate = inf;
end
% Adjust the amount injected by the number of layers used in the model
W(5).val = 5000*stb/(day*numel(layers));
W(5).lims.bhp = 10000*psia;
W(5).lims.rate = inf;
W(5).type = 'rate';


sys = initADISystem(deck, G, rock, fluid);
sys.nonlinear.cprBlockInvert = true;
sys.nonlinear.cpr = true;
sys.nonlinear.cprEllipticSolver = @(A, b) agmg(A, b);
% sys.nonlinear.cprEllipticSolver = @(A, b) mldivide(A, b);
sys.nonlinear.maxIterations = 300;
mrstVerbose on
%%
endtime = 1*day;
dt0 =  .1*day/numel(layers);
%% Solve schedule
dt = dt0;
T = 0;
states = [];
state = state0;
dt_history = [];
ithist = [];
timer = tic();
its = 0;
while T < endtime
    [dt, dt_history] = simpleStepSelector(dt_history, dt, its, 'targetIts', 10);
    [state, its, convergence] = solvefiADI(state, dt, W, G, sys);
    T = T + dt;
    states = [states; state];
    ithist = [ithist, its];
    fprintf([formatTimeRange(T), '\n']);
end
time_cpr = toc(timer);
%%
clf;
plotToolbar(G, states)
