function G = computeEffectiveTrans(G)
% Computes the effective transmissibility between the fracture and matrix
% by harmonically averaging weighted permeablities and multiplying by CI.

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

if isfield(G.rock,'poro'), pv = poreVolume(G,G.rock);
else pv = G.cells.volumes; end
w1 = pv(G.nnc.cells(:,1))./G.rock.perm(G.nnc.cells(:,1));
w2 = pv(G.nnc.cells(:,2))./G.rock.perm(G.nnc.cells(:,2));
wt = pv(G.nnc.cells(:,1))+pv(G.nnc.cells(:,2));
% No weighting by PV:
% w1 = 1./G.rock.perm(G.nnc.cells(:,1));
% w2 = 1./G.rock.perm(G.nnc.cells(:,2));
% wt = 1;
G.nnc.T = G.nnc.CI.*(wt./(w1+w2));

return