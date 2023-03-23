function [Gl,T] = makeNNCextruded(G,Gl,F,fracture,flayers)
% makeNNCextruded identifies fracture-matrix and fracture-fracture
% connections and computes a single transmissibility for each connection. 
%
% SYNOPSIS:
%   [Gl,T] = makeNNCextruded(G,Gl,F,fracture,flayers)
%
% DESCRIPTION:
%     defineNNCandTrans first recognizes all fracture-matrix connections in
%     2D and stores them as non-neighboring connections (NNC's) following
%     which these connections are added to the layerd grid. Then, the
%     global grid 'G' containing both fracture and matrix grid cells (in 3D) is
%     assembled. Next, fracture-fracture connections (at
%     intersections) are defined as NNC's and added to the global grid. The
%     script also computes transmissibilities for these NNC's using rock
%     properties defined in G. This function calls frac_matrix_nnc,
%     frac_frac_nnc, assembleGlobalGrid and computeTrans internally.
%
% REQUIRED PARAMETERS:
%
%   G           - Grid data structure containing G.FracGrid (see
%                 FracTensorGrid2D) with corresponding rock properties and
%                 G.cells.fracture (see markcells2D and CIcalculator2D)
%
%   Gl          - Extruded matrix and fracture grids. See makeLayers
%
%   F, fracture - Output from gridFracture2D.
%
%   flayers     - Indices of extruded layers in which fractures are present.
%
% RETURNS:
%   Gl - Global extruded grid structure with matrix and fracture assemblies
%        combined into a single grid. NNC's and other relevant information
%        is defined in G.nnc.
%
%   T  - Face transmissibilities.
%
% SEE ALSO:
%   assembleGlobalGrid, frac_matrix_nc, frac_frac_nnc_extruded, makeLayers 

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

nfl = numel(flayers);

% frac-matrix connections
G = frac_matrix_nnc(G,F,fracture);
% Remove unnecessary connections
if any(G.nnc.CI==0)
    G.nnc.cells = G.nnc.cells(G.nnc.CI~=0,:);
    G.nnc.type(G.nnc.CI==0,:) = [];
    G.nnc.CI = G.nnc.CI(G.nnc.CI~=0);
end
nm = G.cells.num;
findex2D = []; findex3D = [];
for i = 1:numel(F)
    findex2D = [findex2D,F(i).cells.start]; %#ok
    findex3D = [findex3D, ...
        Gl.FracGrid.(['Frac',num2str(i)]).cells.start]; %#ok
end
findex2D(end+1) = findex2D(end) + F(i).cells.num;
findex3D(end+1) = findex3D(end) + Gl.FracGrid.(['Frac',num2str(i)]).cells.num;
Gl.nnc.cells = []; Gl.nnc.CI = repmat(G.nnc.CI,nfl,1);
Gl.nnc.type = repmat(G.nnc.type,nfl,1);
for i = 1:size(G.nnc.cells,1)
    mat = G.nnc.cells(i,1);
    frac = G.nnc.cells(i,2);
    plane = find(frac>=findex2D(1:end-1) & frac<findex2D(2:end));
    ind2D = frac - findex2D(plane);
    ind3D = findex3D(plane) + ind2D;
    temp = [mat + (flayers'-1)*nm, ind3D + ((1:nfl)'-1)*F(plane).cells.num];
    Gl.nnc.cells = [Gl.nnc.cells;temp];
end
    
Gl = assembleGlobalGrid(Gl); 
Gl.fractureAperture = fracture.aperture;
Gl = computeEffectiveTrans(Gl);

Gl = frac_frac_nnc_extruded(G, Gl, F, fracture, flayers);
Gl = rmfield(Gl,{'numLayers','layerSize','type'});
% compute transmissibilities
T = computeTrans(Gl, Gl.rock);
%-------------------------------------------------------------------------%
% computeTrans returns 2 half-transmissibilities for each internal face and
% one transmissibility for each external face. Below we compute one
% transmissibility per face.
%-------------------------------------------------------------------------%
cf = Gl.cells.faces(:,1);
nf = Gl.faces.num;
T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
T = [T;Gl.nnc.T];

Gl.faces.tag = zeros(Gl.faces.num,1);
x = ismember(Gl.faces.neighbors,Gl.nnc.cells,'rows');
Gl.faces.tag(x) = 1; % Tag NNC faces. 
return

function Gl = computeEffectiveTrans(Gl)
if isfield(Gl.rock,'poro'), pv = poreVolume(Gl,Gl.rock);
else pv = Gl.cells.volumes; end
w1 = pv(Gl.nnc.cells(:,1))./Gl.rock.perm(Gl.nnc.cells(:,1));
w2 = pv(Gl.nnc.cells(:,2))./Gl.rock.perm(Gl.nnc.cells(:,2));
wt = pv(Gl.nnc.cells(:,1))+pv(Gl.nnc.cells(:,2));
% No weighting by PV:
% w1 = 1./Gl.rock.perm(Gl.nnc.cells(:,1)); 
% w2 = 1./Gl.rock.perm(Gl.nnc.cells(:,2));
% wt = 1;
Gl.nnc.T = Gl.nnc.CI.*(wt./(w1+w2));
return