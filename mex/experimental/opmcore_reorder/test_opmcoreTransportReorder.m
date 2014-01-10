
% Set up model
G         = computeGeometry(cartGrid([100,100]));
rock.perm = ones(G.cells.num, 1);
rock.poro = rock.perm;
fluid     = initSimpleFluid('mu', [1, 1], 'rho', [1, 1], 'n', [1,1]);
state     = initState(G, [], 0, [0,1]);
src       = addSource([], [1, G.cells.num], [1, -1]);

% Solve pressure once
T         = computeTrans(G, rock);
state     = incompTPFA(state, G, T, fluid, 'src', src);


% Setup for reorder code
source = zeros(G.cells.num, 1);
source(src.cell)=src.rate;
s = state.s(:,1);
T = sum(poreVolume(G, rock)); % Inject one pore volume
n = 100;
for i=1:n,
    state = opmcoreTransportReorder(state, G, rock, source, T/n);
    surf(reshape(state.s(:,1), G.cartDims));
    view(115, 35);
    drawnow;
end