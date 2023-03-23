function G = assembleGlobalGrid(Gm)
% assembleGlobalGrid combines matrix and fracture grids into 1 global grid.
%
% SYNOPSIS:
%   G = assembleFracGrid(Gm)
%
% REQUIRED PARAMETERS:
%   G  - Matrix grid data structure (passed through computeGeometry)
%        containing G.FracGrid as returned by FracTensorGrid2D.
%
% RETURNS:
%   G - Grid structure containing both matrix and fracture cells and
%       information about fracture-matrix NNC connections
%
% NOTE:
%   This function assumes that each fracture grid is represented as a
%   cartesian grid. Underlying matrix grid structure as in 'Gm' can be
%   cartesian or unstructured.
%
% SEE ALSO:
%   FracTensorGrid2D

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

if isfield(Gm,'Matrix')
    G = Gm;
    return
end
G = Gm;
Gm = rmfield(Gm,{'FracGrid','nnc'});
G.Matrix = Gm;
name = 'Frac';
for i = 1:numel(fieldnames(G.FracGrid))
    Gf = G.FracGrid.([name,num2str(i)]);
    cell_start = numel(G.cells.volumes); %max(G.cells.indexMap); %Gf.cells.start-1;
    face_start = numel(G.faces.areas); %G.faces.tag %Gf.faces.start-1;
    node_start = size(G.nodes.coords,1); %Gf.nodes.start-1;
    if isfield(G.cells,'indexMap')
        G.cells.indexMap = [G.cells.indexMap;Gf.cells.indexMap+cell_start];
    end
    fPos = Gf.cells.facePos+max(G.cells.facePos)-1;
    G.cells.facePos = [G.cells.facePos;fPos(2:end)];
    if isfield(G,'cartDims') || (size(G.cells.faces,2)==2 && size(Gf.cells.faces,2)==2)
        G.cells.faces = [G.cells.faces;[Gf.cells.faces(:,1)+face_start, Gf.cells.faces(:,2)]];
    else
        G.cells.faces = [G.cells.faces;Gf.cells.faces(:,1)+face_start];
    end
    G.cells.volumes = [G.cells.volumes;Gf.cells.volumes];
    G.cells.centroids = [G.cells.centroids;Gf.cells.centroids];
    fnbrs = Gf.faces.neighbors;
    for j = 1:numel(fnbrs)
        if fnbrs(j)~=0
            fnbrs(j) = fnbrs(j) + cell_start;
        end 
    end
    G.faces.neighbors = [G.faces.neighbors;fnbrs];
    nPos = Gf.faces.nodePos+max(G.faces.nodePos)-1;
    G.faces.nodePos = [G.faces.nodePos;nPos(2:end)];
    if isfield(G.faces,'tag') && isfield(Gf.faces,'tag')
        if ~isfield(Gf.faces,'tag')
            Gf.faces.tag = zeros(Gf.faces.num,1);
        end
        G.faces.tag = [G.faces.tag;Gf.faces.tag];
    end
    G.faces.nodes = [G.faces.nodes;Gf.faces.nodes+node_start];
    G.faces.areas = [G.faces.areas;Gf.faces.areas];
    G.faces.normals = [G.faces.normals;Gf.faces.normals];
    G.faces.centroids = [G.faces.centroids;Gf.faces.centroids];
    G.nodes.coords = [G.nodes.coords;Gf.nodes.coords];
    G.rock.perm = [G.rock.perm;Gf.rock.perm];
    if isfield(Gf.rock,'poro') && isnumeric(Gf.rock.poro)
        G.rock.poro = [G.rock.poro;Gf.rock.poro];
    end
end
G.cells.num = numel(G.cells.volumes); % numel(indexMap)
G.faces.num = max(G.cells.faces(:));
G.nodes.num = size(G.nodes.coords,1);
return
    
