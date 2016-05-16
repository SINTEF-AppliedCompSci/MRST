function [G,T] = defineNNCandTrans(G,F,fracture)
% defineNNCandTrans identifies fracture-matrix and fracture-fracture
% connections and computes once transmissibility for each connections.
%
% SYNOPSIS:
%   [G,T] = defineNNCandTrans(G,F,fracture)
%
% DESCRIPTION:
%     defineNNCandTrans first recognizes all fracture-matrix connections
%     and stores them as non-neighboring connections (NNC's). Then, the
%     global grid 'G' containing both fracture and matrix grid cells is
%     assembled. Following that, fracture-fracture connections (at
%     intersections) are defined as NNC's and added to the global grid. The
%     script also computes transmissibilities for these NNC's using rock
%     properties defined in G. This function calls frac_matrix_nnc,
%     frac_frac_nnc, assembleGlobalGrid and computeTrans internally.
%
% REQUIRED PARAMETERS:
%
%   G           - Grid data structure containing G.FracGrid (see
%                 FracTensorGrid2D) with corresponding rock properties and
%                 G.cells.fracture (see markcells2D and CIcalculator2D).
%
%   F, fracture - Output from gridFracture2D.
%
% RETURNS:
%   G - Global grid structure (individual matrix and fracture grids
%       assembled into 1 grid) with NNC's and relevant information defined
%       in G.nnc.
%
%   T - Face transmissibilities.
%
% SEE ALSO:
%   gridFracture2D, assembleGlobalGrid, frac_matrix_nc, frac_frac_nnc 

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
% Remove unnecessary connections i.e. those with very low CI
if any(G.nnc.CI<=eps*100)
    G.nnc.cells = G.nnc.cells(G.nnc.CI>eps*100,:);
    G.nnc.type = G.nnc.type(G.nnc.CI>eps*100,:);
    G.nnc.CI = G.nnc.CI(G.nnc.CI>eps*100);
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
% computeTrans returns 2 half-transmissibilities for each internal face and
% one transmissibility for each external face. Below we compute one
% transmissibility per face.
%-------------------------------------------------------------------------%
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
T = [T;G.nnc.T];
return