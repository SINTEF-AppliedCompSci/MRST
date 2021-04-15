function H = removeNodes(G, varargin)
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


  warning('Unsafe for faulted grids');

  % For instance, consider
  % G=processGRDECL(readGRDECL('fourz.grdecl'));
  % H=removeNodes(, 'tolerance', 0.3);
  % plotGrid(H, 2); view(3);

  opt = struct('tolerance', 0);
  opt = merge_options(opt, varargin{:});


  e = sort(faceNodes2FaceEdges(G), 2);
  c = G.nodes.coords(e(:,2), :) - G.nodes.coords(e(:,1),:);
  L = sqrt(sum(c.*c, 2));

  %remove small edges by merging nodes.
  remove = L < opt.tolerance;
  if ~any(remove), H = G;return; end

  [map, remove_nodes] = constructNodeMap(G, e, remove);
  G.faces.nodes = map(G.faces.nodes);

  % remove unused nodes
  G.nodes.coords(remove_nodes,:) = [];
  G.nodes.num = size(G.nodes.coords, 1);

  % remove collapsed edges from grid
  [G.faces.nodes, G.faces.nodePos] = removeFromPackedData(G.faces.num, ...
     G.faces.nodePos, G.faces.nodes, remove);

  % remove collapsed faces
  faces    = diff(G.faces.nodePos) < 3;
  G        = gridutils(G, 'remove faces', faces);

  % remove collapsed cells.

  H = G;
  H.type = [H.type, mfilename];
end

function [map, remove] = constructNodeMap(G, e, remove)
   % Hint:
   %
   % r defines an equivalence relation between nodes: a row in r
   % corresponds to a pair of nodes that should be collapsed to one. The
   % following routines determine the new numbering of nodes.
   %
   %
   r  = double (e(remove, :));
   n  = max(r(:));
   A  = sparse(r(:,1), r(:,2), 1, n, n); A=A' + A + speye(n);
   [a,b,c,d] = dmperm(A);

   map    = (1:G.nodes.num)';
   i = rldecode(c(1:end-1), diff(c), 2) .';
   map(a) = a(i);
   renumber=cumsum(accumarray(map, 1)>0);
   map = renumber(map);
   remove = a(mcolon(c(1:end-1)+1, c(2:end)-1)) .';

end

function [r, p] = reverseLookup(n, pos, data, v)
   i    = false(max(data(:,1)), 1);
   i(v) = true;
   no   = rldecode(1:n, diff(pos), 2).';
   p    = i(data(:,1));
   r    = no(p);
end

function [data, pos] = removeFromPackedData(n, pos, data, remove)
   % E.g. [cellFaces, facePos]= removeFromPackedData(cells.num, facePos,...
   %                               cellFaces, faces_to_remove)
   [r, p]  =  reverseLookup(n, pos, data, remove);
   data(p,:) = [];
    pos     = int32(cumsum([1; double(diff(pos))-accumarray(r,1,[n,1])]));
end

function e = faceNodes2FaceEdges(G)
  numnodes = double(rldecode(diff(G.faces.nodePos), ...
                             diff(G.faces.ndoePos)));

  pos      = cumsum([1; double(diff(G.faces.nodePos))]);
  pos      = rldecode(pos(1:end-1), G.faces.numNodes);

  this     = mcolon(ones([size(G.faces.nodePos,1), 1]), ...
                    double(diff(G.faces.nodePos)))' - 1;

  e(:,1) = G.faces.nodes(pos + mod(this ,    numnodes));
  e(:,2) = G.faces.nodes(pos + mod(this + 1, numnodes));
end
