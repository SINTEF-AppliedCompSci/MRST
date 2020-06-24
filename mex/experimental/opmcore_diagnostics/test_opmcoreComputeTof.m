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

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
