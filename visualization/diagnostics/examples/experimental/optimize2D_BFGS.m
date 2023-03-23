mrstModule add ad-fi diagnostics ad-props spe10 incomp optimization
%clear functions

useBasis = false;
useAGMG  = false;

[G, W, rock] = getSPE10setup(75);
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

% set up linear solvers with timing
clear linsolveWithTimings
ls_agmg = @(A,x)linsolveWithTimings(A,x,@agmg);
ls      = @linsolveWithTimings;

if useBasis
    ls_bas = ls;
    if useAGMG, ls_bas = ls_agmg; end
    basis = computeDefaultBasis([], G, state0, s, W, fluid_ad, pv, T, 'linsolve', ls_bas);
else
    basis = [];
end
% get timing:
t_basis = linsolveWithTimings;

ls_pres = ls;
if useAGMG&&~useBasis
    ls_pres = ls_agmg;
end
% set up objective evalueation
f = @(u)evalObjectiveDiagnostics(u, objective, state0, s, G, fluid_ad, pv, T, W, scaling, ...
        'targets', targets, 'linSolve', ls_pres, 'linsolveTOF', ls, 'msbasis', basis);

% set initial controls and run optimization
uInit = well2control(W, 'scaling', scaling, 'targets', targets);
[v, u, hst] = unitBoxBFGS(uInit, f);

% get timing
t_lin = linsolveWithTimings;
t_lin = reshape(t_lin, [6, numel(t_lin)/6]).';

% display optimal well-controls
Wopt = control2well(u, W, 'scaling', scaling, 'targets', targets);
for k = 1:numel(targets)
    wnr = targets(k);
    fprintf([W(wnr).name, ',  initial bhp: %5.1f ,  optimal bhp: %5.1f\n'], W(wnr).val/barsa, Wopt(wnr).val/barsa);
end
% display timing info
fprintf('Timing basis computations: %10.2e s\n', sum(t_basis));
fprintf('Mean linear solver times: \n');
cases = {'Pressure forward            :',...
         'Forward TOF                 :',...
         'Backward TOF                :',...
         'Forward TOF - adjoint       :',...
         'Backward TOF - adjoint      :',...
         'Pressure -adjoint           :'};
mt_lin = mean(t_lin);
for cnr = [1 2 3 5 4 6];
    fprintf([cases{cnr}, '%10.2e', ' s\n'], mt_lin(cnr));
end
return

