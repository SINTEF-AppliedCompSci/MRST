function [G,T] = defineNNC_computeTrans(G,F,fracture)
% defineNNC_computeTrans uses a global grid 'G' with fractures to define
% fracture-matrix connections and fracture-fracture connections at
% intersections as non-neighboring connections (NNC's). The script also
% computes transmissibilities for these NNC's using rock properties defined
% in G. This function calls frac_matrix_nnc, frac_frac_nnc, 
% assembleGlobalGrid  and computeTrans internally.

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


% frac-matrix connections
G = frac_matrix_nnc(G,F,fracture);
% Remove unnecessary connections
if any(G.nnc.CI==0)
    G.nnc.cells = G.nnc.cells(G.nnc.CI~=0,:);
    G.nnc.type(G.nnc.CI==0,:) = [];
    G.nnc.CI = G.nnc.CI(G.nnc.CI~=0);
end

% global grid
G = assembleGlobalGrid(G); 
G.fractureAperture = fracture.aperture;
G = computeEffectiveTrans(G);

% star-delta at fracture intersections
G = frac_frac_nnc(G,F,fracture);

% compute transmissibilities
T = computeTrans(G, G.rock);
%-------------------------------------------------------------------------%
% computeTrans returns 2 transmissibilities for each internal face and one
% transmissibility fo each external face. Size of transmissibility array is
% same as G.cells.faces if opt.usetrans = false 
%-------------------------------------------------------------------------%
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
T = [T;G.nnc.T];
return

function G = computeEffectiveTrans(G)
w1 = G.cells.volumes(G.nnc.cells(:,1))./G.rock.perm(G.nnc.cells(:,1));
w2 = G.cells.volumes(G.nnc.cells(:,2))./G.rock.perm(G.nnc.cells(:,2));
wt = G.cells.volumes(G.nnc.cells(:,1))+G.cells.volumes(G.nnc.cells(:,2));
G.nnc.T = G.nnc.CI.*(wt./(w1+w2)); clear wt w1 w2
return