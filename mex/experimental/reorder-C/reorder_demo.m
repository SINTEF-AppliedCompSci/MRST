%% Some pressure solvers may include the effect of regions, by the hack
%% "call fluid.relperm with as many saturations as there are cells in G".

G = cartGrid([100,100], [10,10]);
G = computeGeometry(G);

rock.perm  = logNormLayers(G.cartDims, 'a', 3);
rock.poro  = ones(G.cells.num, 1);

%% Make special fluid:
satnum = ones(G.cartDims);
satnum(1:end/2, end/2+1:end)     = 2;
satnum(end/2+1:end, 1:end/2)     = 3;
satnum(end/2+1:end, end/2+1:end) = 4;
satnum = satnum(:);

fluid      = initCoreyFluidROC('mu',  [1,1], ...
                               'rho', [0,0], ...
                               'sr',  [0,0;  ...
                                       0,0;  ...
                                       0,0;  ...
                                       0,0], ...
                               'n',   [1,1;  ...
                                       2,2;  ...
                                       2,3;  ...
                                       2,4], ...
                               'kwm', [1,1;  ...
                                       1,1;  ...
                                       1,1;  ...
                                       1,1], ...
                               'reg', satnum);

%% Compute pressure...
disp('Computing transmissibilities...')
T          = computeTrans(G, rock);
S          = computeMimeticIP(G, rock);
state      = initState(G, [], 0);

src        = addSource([], [1, G.cells.num], [1, -1], 'sat', [1,1]);
disp('Computing pressure...')
%state      = incompTPFA(state, G, T, fluid, 'src', src);
state      = solveIncompFlow(state, G, S, fluid, 'src', src);

%% Compute transport
dt = repmat(0.1, [10, 1]);

disp('Computing transport...')
cla,
for i=1:100,
   state = implicitupwind(state, G, dt, rock, fluid, 'src', src);
   %plotCellData(G, state.s(:,1)); axis equal tight
   surf(reshape(state.s(:,1), G.cartDims));shading flat;
   set(gca, 'dataaspectr', [1,1,0.02]);view(130, 10)
   drawnow;
end
