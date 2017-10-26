mrstModule add msrsb incomp coarsegrid spe10
%%
simple = true;
if simple
    % Create n by n grid
    n = 500;
    G = cartGrid([n, n, 1]);
    G = computeGeometry(G);
    rock = makeRock(G, 1, 1);
    W = [];
    W = addWell(W, G, rock, 1, 'type', 'rate', 'val', 1);
    W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 0);
else
    % Take part of SPE10
    layers = 1:35;
    [G, W, rock] = SPE10_setup(layers);
end

state = initResSol(G, 0);
state = incompTPFA(state, G, computeTrans(G, rock), initSingleFluid('rho', 0, 'mu', 1), 'Wells', W, 'matrixOutput', true, 'LinSolve', @(A, b) 0*b);
A = state.A;
b = state.rhs;

%%
try
    mrstModule add agmg
    tic();
    x_agmg = agmg(A, b);
    t_agmg = toc();
catch
    x_agmg = nan*b;
    t_agmg = nan;
end

%% Call AMG version
tic();
x = callAMGCL(A, b, 'coarsening', 'smoothed_aggregation',...
                    'relaxation', 'spai0', ...
                    'preconditioner', 'amg');
t_amg = toc();
%%

tic();
x = callAMGCL(A, b, 'relaxation', 'ilu0',...
                    'preconditioner', 'relaxation');
t_relax = toc();
%%
bar([t_amg, t_relax, t_agmg])
set(gca, 'XTicklabel', {'AMGCL-amg', 'AMGCL-relax' 'AGMG'});
ylabel('Solution time [s]');
%%
