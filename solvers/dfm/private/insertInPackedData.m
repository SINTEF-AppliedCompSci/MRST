function [data, pos] = insertInPackedData(pos, data, r, c)
%Insert c into row r of packed array (data, pos)
%
% SYNOPSIS:
%   [data, pos] = insertInPackedData(n, pos, data, r, c)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
%   fnodes  - face nodes
%
%   nnodes  - number of nodes for each face
%   neigh   - cell neighbors of each face
%   tags    - cellFaces tags for each half-face associated with face.
%
% RETURNS:
%   G       - Grid structure where the followinf fields have been modified:
%
%                  G.faces.nodes
%
%                  G.cells.faces
%                  G.cells.facePos
%                  G.cells.numFaces
%
%                  G.faces.neighbors
%                  G.faces.numNodes
%                  G.faces.nodePos
%                  G.faces.tag
%                  G.faces.num
%
% COMMENTS:
%  Cells may not be closed polygons/polyhedra after this call.
%
% SEE ALSO:

% E.g. [cellFaces, facePos]= insertInPackedData(cells.num, facePos,...
   %                               cellFaces, cells, [facenums, tags])

   n   = numel(pos)-1;
   tmp = sortrows([r, c]);
   new = tmp(:,2:end);
   t   = accumarray(r, 1, [n, 1]);
   num = diff(pos);
   pos = int32(cumsum([1;double(diff(pos))+t]));

   r   = unique(r);
   ix  = mcolon(pos(r)+num(r), pos(r+1)-1);
   newData = zeros(size(data, 1)+size(new, 1), size(data, 2));
   i    = (1:size(newData))';
   i(ix)=[];

   newData(i, :) = data;
   newData(ix,1:size(new,2)) = new;
   data = newData;

