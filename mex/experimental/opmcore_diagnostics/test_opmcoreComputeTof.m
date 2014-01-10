% Set up model
mrstModule add spe10

[G, W, rock] = SPE10_setup(25);
rock.poro = max(rock.poro, 1e-4);
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
state = initState(G, W, 0);
S  = computeMimeticIP(G, rock);
state = solveIncompFlow(state, G, S, fluid, 'wells', W);
source = zeros(G.cells.num, 1);
source(vertcat(W.cells)) = vertcat(state.wellSol.flux);

% Run tof code and display results
tof = opmcoreComputeTof(state, G, rock, source);
surf(reshape(log(tof), G.cartDims));
