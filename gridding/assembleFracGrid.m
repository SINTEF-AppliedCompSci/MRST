function Gf = assembleFracGrid(G)
% assembleFracGrid assembles a grid (Gf) containing only fractures using
% the global grid (G) containing both fracture and matrix information.
%
% SYNOPSIS:
%   Gf = assembleFracGrid(G)
%
% REQUIRED PARAMETERS:
%   G  - Grid structure with fractures (see assembleGlobalGrid)
%
% RETURNS:
%   Gf - Grid structure containing all fracture lines as 1 grid. 
%
% NOTE:
%   This function assumes that each fracture grid is represented as a
%   cartesian grid. Underlying matrix grid structure as in 'G' can be
%   cartesian or unstructured. Fracture intersections exist as NNC's
%   (Gf.nnc.cells) in the returned grid Gf. Transmissibility at these
%   intersections can also be found in Gf.nnc. Fracture-matrix connections
%   do not exist in the output grid.
%
% SEE ALSO:
%   assembleGlobalGrid

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

FG = G.FracGrid;
for i = 1:numel(fieldnames(FG))
    if i == 1, Gf = FG.(['Frac',num2str(i)]);
    else
        Gl = FG.(['Frac',num2str(i)]);
        cell_start = max(Gf.cells.indexMap); %Gf.cells.start-1;
        face_start = numel(Gf.faces.areas); %Gf.faces.start-1;
        node_start = size(Gf.nodes.coords,1); %Gf.nodes.start-1;
        Gf.cells.indexMap = [Gf.cells.indexMap;Gl.cells.indexMap+cell_start];
        fPos = Gl.cells.facePos+max(Gf.cells.facePos)-1;
        Gf.cells.facePos = [Gf.cells.facePos;fPos(2:end)];
        
        if isfield(Gf,'cartDims') || size(G.cells.faces,2)==2
            Gf.cells.faces = [Gf.cells.faces;[Gl.cells.faces(:,1)+face_start, Gl.cells.faces(:,2)]];
        else
            Gf.cells.faces = [Gf.cells.faces;Gl.cells.faces(:,1)+face_start];
        end
        
        Gf.cells.volumes = [Gf.cells.volumes;Gl.cells.volumes];
        Gf.cells.centroids = [Gf.cells.centroids;Gl.cells.centroids];
        fnbrs = Gl.faces.neighbors;
        for j = 1:numel(fnbrs)
            if fnbrs(j)~=0
                fnbrs(j) = fnbrs(j) + cell_start;
            end
        end
        Gf.faces.neighbors = [Gf.faces.neighbors;fnbrs];
        nPos = Gl.faces.nodePos+max(Gf.faces.nodePos)-1;
        Gf.faces.nodePos = [Gf.faces.nodePos;nPos(2:end)];
        if isfield(Gf.faces,'tag')
            Gf.faces.tag = [Gf.faces.tag;Gl.faces.tag];
        end
        Gf.faces.nodes = [Gf.faces.nodes;Gl.faces.nodes+node_start];
        Gf.faces.areas = [Gf.faces.areas;Gl.faces.areas];
        Gf.faces.normals = [Gf.faces.normals;Gl.faces.normals];
        Gf.faces.centroids = [Gf.faces.centroids;Gl.faces.centroids];
        Gf.nodes.coords = [Gf.nodes.coords;Gl.nodes.coords];
        if isfield(Gf,'rock') 
            Gf.rock.perm = [Gf.rock.perm;Gl.rock.perm];
            if isfield(Gf.rock,'poro') && isnumeric(Gf.rock.poro)
                Gf.rock.poro = [Gf.rock.poro;Gl.rock.poro];
            end
        end
    end
    Gf.cells.num = max(Gf.cells.indexMap);
    Gf.faces.num = max(Gf.cells.faces(:));
    Gf.nodes.num = size(Gf.nodes.coords,1);
    Gf.cartDims = [Gf.cells.num 1];
    if isfield(G.nnc,'type')
        Gf.nnc.cells = G.nnc.cells(strcmp(G.nnc.type,'star-delta'),:)-G.Matrix.cells.num;
        Gf.nnc.T = G.nnc.T(strcmp(G.nnc.type,'star-delta'));
    end
end