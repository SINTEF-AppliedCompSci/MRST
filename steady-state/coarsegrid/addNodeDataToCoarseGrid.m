function CG = addNodeDataToCoarseGrid(CG)
% Add node data to upscaled grid CG
% 
% Assumptions:
% - Each face has four nodes.
% - Each face is a rectangle
% - Each face is normal to one of the main directions
% 
% NOTE: NOT SURE ABOUT THE ORDERING OF THE NODES WITHIN EACH FACE. THIS MAY
% BE A PROBLEM IF GRID IS PASSED ON TO E.G. computeGeometry, WHICH ASSUMES
% A SPEECIAL ORDERING OF THE NODES.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

G = CG.parent;

% Allocate space for the nodes
faceNodes = nan(CG.faces.num, 4);

% Loop over coarse scale faces
for i = 1:CG.faces.num
   
   % Fine scale faces within the coarse scale face
   faces = CG.faces.fconn(CG.faces.connPos(i) : ...
                          CG.faces.connPos(i + 1) - 1);
   
   % Allocate space for fine scale nodes
   nodes = nan(numel(faces)*4,1);
   
   % Loop over fine scale faces
   for j = 1:numel(faces)
      
      % Face index
      f = faces(j);
      
      % Get fine scale nodes for this face
      fnodes = G.faces.nodes( G.faces.nodePos(f) : ...
         G.faces.nodePos(f+1)-1 );
      assert(numel(fnodes) == 4, ...
         'Assumption that each face has four nodes does not hold.');
      nodes((j-1)*4 + 1 : j*4) = fnodes;
      
   end
   
   % Get unique nodes and their coords
   nodes  = unique(nodes);
   coords = G.nodes.coords(nodes, :);
   
   % Direction of normal to face
   [~, ndir] = max(abs(CG.faces.normals(i,:)), [], 2);
   
   % Extract only the two 'main' directions
   dims = [1:ndir-1, ndir+1:3];
   coordsm = coords(:, dims);
   
   tol = 1e-6;
   minmax = [min(coordsm, [], 1); max(coordsm, [], 1)];
   ext = @(m1, m2) find(abs(coordsm(:,1) - minmax(m1, 1)) < tol & ...
                        abs(coordsm(:,2) - minmax(m2, 2)) < tol);
   corners = [ext(1,1); ext(2,1); ext(2,2); ext(1,2)]';
   
   % TODO: NOT SURE ABOUT THE ORDER. SEE ALSO tiltedGrid (174).
   
   % Extract corners info
   faceNodes(i, :) = nodes(corners)';
   
end

faceNodes       = reshape(faceNodes', [], 1);
faceNodesUnique = unique(faceNodes);
faceNodesLocal  = nan(numel(faceNodes), 1);
for i = 1:numel(faceNodes)
   faceNodesLocal(i) = find(faceNodes(i) == faceNodesUnique);
end

% Add data to grid structure
CG.faces.nodePos = (1:4:(4*CG.faces.num)+1)';
CG.faces.nodes   = faceNodesLocal;
CG.nodes.num     = numel(faceNodesUnique);
CG.nodes.coords  = G.nodes.coords(faceNodesUnique, :);
CG.nodes.global  = faceNodesUnique;

end

