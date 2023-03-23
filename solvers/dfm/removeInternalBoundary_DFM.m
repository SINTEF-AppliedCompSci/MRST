function G = removeInternalBoundary_DFM(G, N)
%Remove internal boundary in grid by merging faces in face list N
%
% This function is modified to account for hybrid cells.
%
% SYNOPSIS:
%   G = removeInternalBoundary_DFM(G, N)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
%   N       - An n x 2 array of face numbers.  Each pair in the array
%             will be merged to a single face in G.  The connectivity if
%             the new grid is updated accordingly.  The geometric adjacency
%             of faces is not checked.
%
% RETURNS:
%   G       - Modified grid structure.
%
%
% COMMENTS:
%
%  What if nodes in f1 are permuted compared to nodes in f2, either due to
%  sign of face (2D) or due to arbitary starting node (3D)?  For the time
%  being, this code assumes that nodes that appear in faces that are being
%  merged coincide exactly --- no checking is doen on node positions.  If
%  nodes are permuted in one of the faces, the resulting grid will be
%  warped.
%
% SEE ALSO:
%  `makeInteralBoundary`


%{
Copyright 2009, 2010, 2011 SINTEF ICT, Applied Mathematics.

Portions Copyright 2011-2012 University of Bergen.

This file is part of DFM module of The MATLAB Reservoir Simulation Toolbox
(MRST).

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


if isempty(N)
    return;
end

N = N(N(:,1) ~= N(:,2), :);
if isempty(N),
    warning('No internal boundaries removed');%#ok
    return;
end
   % Save faces on each side of boundary
   f1 = N(:,1);
   f2 = N(:,2);

   % Copy face info  (sequence of stuff is uncertain: check and fix!)
   [fnodes1, nnodes1, neigh1, tags1] = copyFaces(G, f1);
   [fnodes2, nnodes2, neigh2, tags2] = copyFaces(G, f2); %#ok

   % Replace nodes in fnodes2 by nodes in fnodes1
    G = replaceNodes(G, fnodes1, fnodes2);

   % Should remove unused nodes.
    %used_nodes = accumarray(G.faces.nodes, 1, [G.nodes.num,1]);

   % Remove all faces
   [G,map] = removeFaces(G, N(:));



   % Merge face data, preserve sign of f1.
   i        = neigh1' == 0;
   if any(sum(i)==2)
       i(1,sum(i)==2)=0;
   end
   neigh    = neigh1';
   neigh(i) = sum(neigh2, 2);
   neigh    = neigh';




   %tags     = tags1';
   %tags(i)  = sum(tags2, 2);
   %tags     = tags';

   % Add merged faces.  Make sure this is f1, since we  removed node
   % numbers associated with f2.
   G = addFaces(G, fnodes1, nnodes1, neigh);

   % Also update face properties;
   % Use f1 since thats the nodes we keep
   if isfield(G.faces,'areas')
       areas = G.faces.areas(N(:,1));
       %areas = mean(G.faces.areas(N),2);
       G.faces.areas = [G.faces.areas(map>0); areas];
   end
   if isfield(G.faces,'normals')
       normals = G.faces.normals(N(:,1),:);
       G.faces.normals = [G.faces.normals(map>0,:);normals];
   end
   if isfield(G.faces,'centroids')

       centroids = G.faces.centroids(N(:,1),:);
       G.faces.centroids = [G.faces.centroids(map>0,:); centroids];
   end
   if isfield(G.faces,'tags')
       tags = G.faces.tags(N(:,1));
       G.faces.tags =[G.faces.tags(map>0); tags];
   end

   if isfield(G.faces, 'hybrid')
       hybrid = G.faces.hybrid(N(:,1));
       G.faces.hybrid = [G.faces.hybrid(map>0); hybrid];
   end

    if isfield(G.faces, 'level')
       level = G.faces.level(N(:,1));
       G.faces.level = [G.faces.level(map>0); level];
    end
     if isfield(G.faces, 'coarse')
        coarse = G.faces.coarse(N(:,1));
        G.faces.coarse = [G.faces.coarse(map>0); coarse];
     end

    if isfield(G,'hybridNeighbors')
        faces = G.hybridNeighbors.faces;
        if ~isempty(faces)
            G.hybridNeighbors.faces  = map(faces);
        end
    end

    if isfield(G.cells, 'tags2')

        G.cells.tags2(G.cells.tags2>0) = map(G.cells.tags2(G.cells.tags2>0));
    end
   G.type = [G.type, { mfilename }];
end

function G = replaceNodes(G, replacement, remove)
   replacement = flipNodes(G,replacement,remove);
   map           = (1:G.nodes.num)';
   map(remove)   = replacement;
   G.faces.nodes = map(G.faces.nodes);
   G = removeNodes(G,remove);
end

function replacement = flipNodes(G,replacement,remove)

    % make sure that the removed nodes has the same ordering
    remove = reshape(remove,2,[])';
    replacement = reshape(replacement,2,[])';

    %flip if the vector from node 1 to node 2 do not point in the same direction.
    doflip = dot(G.nodes.coords(remove(:,1),:)-G.nodes.coords(remove(:,2),:),G.nodes.coords(replacement(:,1),:)-G.nodes.coords(replacement(:,2),:),2)<0;
    replacement(doflip,:) = fliplr(replacement(doflip,:));
    replacement = reshape(replacement',[],1);
end
function G = removeNodes(G,remove)
   map  = (1:G.nodes.num)';
   map2  = setdiff(map,remove);
   G.nodes.coords = G.nodes.coords(map2,:);
   G.nodes.tags = G.nodes.tags(map2,:);
   G.nodes.level = G.nodes.level(map2,:);
   G.nodes.num = G.nodes.num - numel(unique(remove));
   map(map2) = (1:G.nodes.num)';
   G.faces.nodes = map(G.faces.nodes);


end
