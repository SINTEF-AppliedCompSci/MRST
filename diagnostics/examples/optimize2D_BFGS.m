mrstModule add ad-fi diagnostics ad-props spe10

[G, W, rock] = SPE10_setup(1);
rock.poro = max(rock.poro, 1e-4);


T = computeTrans(G, rock);
pv = poreVolume(G, rock);

% Reservoir state
state0 = initResSol(G, 300*barsa, [1 0 0]);

% Reservoir fluid and ad-system
mu = [1 1 1];
n = [1 1 1];

fluid_ad = initSimpleADIFluid('mu', mu, 'n', n);

s = initADISystem({'Oil', 'Water'}, G, rock, fluid_ad);

objective = getObjectiveDiagnostics(G, rock, 'minlorenz');
% only control well 1:4:
targets = (1:4)';

% set limits/scaling (don't scale objective, already ~O(1)),
% producers alowed to operate between 100 and 300 bar
scaling.boxLims = repmat([100 300], numel(targets), 1)*barsa;
scaling.obj     = 1;

% set up objective evalueation
f = @(u)evalObjectiveDiagnostics(u, objective, state0, s, G, fluid_ad, pv, T, W, scaling, 'targets', targets);

% set initial controls and run optimization
uInit = well2control(W, 'scaling', scaling, 'targets', targets);
[v, u, hst] = unitBoxBFGS(uInit, f);

% display optimal well-controls
Wopt = control2well(u, W, 'scaling', scaling, 'targets', targets);
for k = 1:numel(targets)
    wnr = targets(k);
    fprintf([W(wnr).name, ',  initial bhp: %5.1f ,  optimal bhp: %5.1f\n'], W(wnr).val/barsa, Wopt(wnr).val/barsa);
end
return

