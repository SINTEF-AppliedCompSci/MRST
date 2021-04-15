%run ../../startup
G = computeGeometry(cartGrid([200e3, 1], [1, 1]));
rock.perm = ones(G.cells.num, 1);

BI = mex_ip_simple(G, rock);

connPos = G.cells.facePos;
conns   = G.cells.faces(:,1);

[S, r, F, L, q] = mex_schur_comp_symm(BI, connPos, conns, []);

nconn  = diff(connPos);
[i, j] = blockDiagIndex(nconn, nconn);

SS = sparse(double(conns(i)), double(conns(j)), S);
R  = accumarray(conns, r);

SS(1) = SS(1) * 2;
R([1, G.cells.num+1]) = [1, -1];

x = SS \ R;
%%
[v, p] = mex_compute_press_flux(BI, x, connPos, conns, F, L, q);

%%
plotCellData(G, p);

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
