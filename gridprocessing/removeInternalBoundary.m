function [G, f] = removeInternalBoundary(G, N)
%Remove internal boundary in grid by merging faces in face list N
%
% SYNOPSIS:
%   G = removeInternalBoundary(G, N)
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
%   f       - New face numbers for the faces that have been merged.
%
%
% NOTE:
%
%  What if nodes in f1 are permuted compared to nodes in f2, either due to
%  sign of face (2D) or due to arbitary starting node (3D)?  For the time
%  being, this code assumes that nodes that appear in faces that are being
%  merged coincide exactly --- no checking is done on node positions.  If
%  nodes are permuted in one of the faces, the resulting grid will be
%  warped.
%
% SEE ALSO:
%  `makeInternalBoundary`


%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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


   if isempty(N),
      warning('No internal boundaries removed');%#ok
      f = [];
      return;
   end

   N = uniqueStable(N(~any(isnan(N), 2) & (N(:,1) ~= N(:,2)), :), 'rows');
   % Save faces on each side of boundary
   f1 = N(:,1);
   f2 = N(:,2);

   % Copy face info  (sequence of stuff is uncertain: check and fix!)
   [fnodes1, nnodes1, neigh1, tags1] = copyFaces(G, f1);
   [fnodes2, nnodes2, neigh2, tags2] = copyFaces(G, f2); %#ok

   % Replace nodes in fnodes2 by nodes in fnodes1
   G = replaceNodes(G, fnodes1, fnodes2);

   % Should remove unused nodes.
   % used_nodes = accumarray(G.faces.nodes, 1, [G.nodes.num,1]);

   % Remove all faces
   G = removeFaces(G, N(:));

   % Merge face data, preserve sign of f1.
   i        = neigh1' == 0;
   neigh    = neigh1';
   neigh(i) = sum(neigh2, 2);
   neigh    = neigh';

   if any(tags1(:)),
      tags     = tags1';
      tags(i)  = sum(tags2, 2);
      tags     = tags';
   end
   % Add merged faces.  Make sure this is f1, since we  removed node
   % numbers associated with f2.
   if any(tags1(:)),
      [G, f] = addFaces(G, fnodes1, nnodes1, neigh, tags);
   else
      [G, f] = addFaces(G, fnodes1, nnodes1, neigh);
   end
   G.type   = [G.type, { mfilename }];
end

%--------------------------------------------------------------------------

function G = replaceNodes(G, replacement, remove)
   map           = (1:G.nodes.num)';
   map(remove)   = replacement;
   G.faces.nodes = map(G.faces.nodes);
end
