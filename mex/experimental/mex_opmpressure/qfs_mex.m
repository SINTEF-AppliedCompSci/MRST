cartDims = [500, 500];
physDims = [1000, 1000];
g = computeGeometry(cartGrid(cartDims, physDims));

rock.perm = ones(g.cells.num, 1);

BI = mex_ip_simple(g, rock);

%{
W = addWell([], g, rock, 1);
W = addWell(W , g, rock, g.cells.num);

[BI, connPos, conns] = mex_ip_simple(g, rock, W);

nconn = diff(connPos);

BI([nconn(1)^2, end]) = [ W.WI ];

[S, r, F, L] = mex_schur_comp_symm(BI, connPos, conns);

[i, j] = blockDiagIndex(nconn, nconn);
SS = sparse(double(conns(i)), double(conns(j)), S);
R  = accumarray(conns, r);

R(end-1 : end) = 1000 * [1, -1];

lam = SS \ R;

[flux, press] = mex_compute_press_flux(BI, lam, connPos, conns, F, L);

plotCellData(g, press);
%}

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
