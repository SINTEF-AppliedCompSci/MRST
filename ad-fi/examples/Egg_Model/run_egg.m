require ad-fi deckformat

%addpath ../../wells_ad-fi
fn    = 'BENCH_EGG.DATA';

deck = readEclipseDeck(fn);

% The deck is given in field units, MRST uses metric.
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

% Create a special ADI fluid which can produce differentiated fluid
% properties.
fluid = initDeckADIFluid(deck);

% The case includes gravity
gravity on

%%
%% Approximate initial conds:
pr   = 400*barsa;
rz   = G.cells.centroids(1,3);
dz   = G.cells.centroids(:,3) - rz;
rhoO    = fluid.bO(400*barsa)*fluid.rhoOS;
rhoW    = fluid.bW(400*barsa)*fluid.rhoWS;
rhoMix  = .1*rhoW + .9*rhoO;
p0   = pr + norm(gravity)*rhoMix*dz;


rSol  = initResSol(G, p0, [0.1, .90]);
%%
system = initADISystem(deck, G, rock, fluid, 'cpr', true);
%system.nonlinear.adhocSolve = true;
system.pscale = 1/(100*barsa);
system.nonlinear.cprBlockInvert = false;
system.nonlinear.cprRelTol      = 2e-2;
system.nonlinear.cprEllipticSolver = @mldivide;
%system.nonlinear.relaxation  = false;
% 
schedule = deck.SCHEDULE;
tt = tic;
wellSols = runScheduleADI(rSol, G, rock, system, schedule);
toc(tt)
% plot 

[wrt, ort, grt, bhp] = wellSolToVector(wellSols);
T = convertTo(cumsum(deck.SCHEDULE.step.val), day);
figure(1), hold on
plot(T,ort(:,9:12)*day), plot(T,wrt(:,9:12)*day),legend(W(9:12).name)

%--------------------------------------------------------------------------
% Timings 21.10.2013:
%   Direct solver:          784 sec
%   CPR (DRS and tol 2e-2)  624 sec
%   Ad-hoc solver           500 sec

