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
for i=1:n
    state = opmcoreTransportReorder(state, G, rock, source, T/n);
    surf(reshape(state.s(:,1), G.cartDims));
    view(115, 35);
    drawnow;
end

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
