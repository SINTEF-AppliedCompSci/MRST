%{
G = computeGeometry(cartGrid([30,30,1]));
src = [];
src = addSource(src, 1, 1);
src = addSource(src, G.cells.num, -1);
bc = [];
rock.perm = ones(G.cells.num, 1);
rock.poro = ones(G.cells.num, 1);

W = addWell([], G, rock, 1, 'type', 'bhp', 'val', 200*barsa);

x = initResSol(G, 0, 0);
x = mex_ifsh(x, G, rock, W, bc, src)

plotCellData(G, x.pressure, 'edgec', 'k', 'edgea', .1, 'facea', .625)
view(3), axis tight, grid on, colorbar southoutside
%}
%%{
clear
G = computeGeometry(cartGrid([3, 1]));
rock.perm = ones([G.cells.num, 1]);

%bc = pside([], G, 'left', 1);
bc = [];
bc = pside(bc, G, 'right', 0);
W  = [];
W  = addWell(W, G, rock, 1, 'type', 'bhp', 'val', 1);

x = initResSol(G, 0);

[x, wbhp, wflux] = mex_ifsh(x, G, rock, W, bc, []);
%}

x2 = initState(G, W, 0);

x2 = solveIncompFlow(x2, G, computeMimeticIP(G, rock),   ...
                     initSingleFluid('mu', 1, 'rho', 1), ...
                     'wells', W, 'bc', bc);

fprintf('Cell pressure error        : %12.5e [relative]\n', ...
        norm(x2.pressure - x.pressure, inf) / ...
        norm(x2.pressure             , inf));

fprintf('Face flux error            : %12.5e [relative]\n', ...
        norm(x2.flux - x.flux, inf) / norm(x2.flux, inf));

fprintf('Well perforation flux error: %12.5e [relative]\n', ...
        norm(vertcat(x2.wellSol.flux) - wflux, inf) / ...
        norm(vertcat(x2.wellSol.flux)        , inf));

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
