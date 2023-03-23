function [G, cellmap] = removeShortEdges(G, tol)
%Replace short edges in grid G by a single node.
%
% SYNOPSIS:
%   G = removeShortEdges(G, tol)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
% RETURNS:
%   G       - Grid structure. Nodes joined by an edge with ||edge||<tol
%             are collapsed to a single node at the edge midpoint. Faces
%             and cells that collapse as a consequence are removed.
%   tol     - OPTIONAL
%             The maximum length of edges that are removed. Default: 0.0
%
% NOTE:
%   This is useful to remove very short edges that may appear in
%   PEBI/Voronoi grids and on faults in cornerpoint grids.
%
% SEE ALSO:
%   `pebi`, `processGRDECL`


%{
The MIT License (MIT)

Copyright (c) 2014 Jostein Natvig
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
%}

   if nargin == 1
      tol = 0.0;
   end

   E = double(faceEdges(G.faces.nodes, diff(G.faces.nodePos)));
   %E = unique(sort(E,2), 'rows');
   
   x = G.nodes.coords;
   L = sqrt( sum( (x(E(:,1),:)-x(E(:,2),:)).^2, 2) );
   i = L<tol;
   
   G = removeEdges(G, E(i,:));
   G = removeCollapsedFaces(G);
   [G, cellmap] = removeCollapsedCells(G);
   
end

function G = removeCollapsedFaces(G)
   % Remove repeated nodes in G.faces.nodes
   faceno = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2)';
   tmp     = [faceno, G.faces.nodes];
   [~, ia] = unique(tmp, 'rows');
   keep    = false(size(tmp,1),1); keep(ia) = true;
   tmp     = tmp(keep,:);
   [n,n] = rlencode(tmp(:,1));
   G.faces.nodePos = cumsum([1;n]);
   G.faces.nodes = tmp(:,2);
   
   % remove collapsed faces
   f = find(diff(G.faces.nodePos)<G.griddim);
   if any(f),G = removeFaces(G, f);end
   
end
   
function [G, cellmap] = removeCollapsedCells(G)
    % Find collapsed faces
   c = find(diff(G.cells.facePos)<G.griddim+1);
   if ~any(c)
       cellmap = (1:G.cells.num)';
       return
   end
   % Remove collapsed cells
   [G, cellmap] = removeCells(G, c);
   % Fix topology
   N = findMatchingFaces(G);
   G = removeInternalBoundary(G, N);   
end

function G = removeEdges(G, E)
% If there are repeated rows in E, it does not matter as dmperm is
% unaffected by values.
%
% This code can be made smarter: only analyze nodes in the set {E(:)}.

   % There can be connected components of graph...
    A = sparse(E(:,1),E(:,2),1,G.nodes.num,G.nodes.num);
    A = A+A'+speye(size(A));
   
    [a,b,c,d]=dmperm(A);

    % find groups of nodes that will be merged
    nno = rldecode(1:numel(c)-1, diff(c), 2)';
    map = zeros(G.nodes.num, 1);
    map(a) = nno;

    x   = accumarray(nno, G.nodes.coords(a,1))./accumarray(nno, 1);
    y   = accumarray(nno, G.nodes.coords(a,2))./accumarray(nno, 1);
    if G.griddim == 3
        z = accumarray(nno, G.nodes.coords(a,3))./accumarray(nno, 1);
        G.nodes.coords = [x,y,z];
    else
        G.nodes.coords = [x,y];
    end
    G.faces.nodes  = map(G.faces.nodes);
    G.nodes.num    = numel(x);
    G.nodes.original_ID = nan;
end


function [E, faceno] = faceEdges(fnodes, nnodes)
% PURPOSE
%
% PARAMETERS
%
% RETURNS
%
% EXAMPLE
%
% SEE ALSO
%
% MISSING
%
% Find edges for each face in G.  Each edge is a pair of node numbers
% (referring to unique nodes in G).  For each face, the edges are listed
% such that face-edges are directed, that is E(i,2) == E(i+1,1), and
% oriented clokckwise. (CHECK THIS)
   nodepos = cumsum([1;double(nnodes(:))]);
   E       = [fnodes, fnodes(rot(nodepos, 1))];
   faceno  = rldecode(1:numel(nnodes), nnodes, 2)';
end

% ------------------------------------------------------------------------

function ix= rot(pos, offset)
% Index rotated OFFSET positions modulo section size in packed array. 
%
% PURPOSE
%
% USAGE
%   Useful operation on packed arrays like (G.faces.nodes,G.faces.nodePos);
%   The index produced by ix = rot(G.faces.nodePos, 1) can be used to
%   rotate nodes for each face by indexing: G.faces.nodes(ix).   
%
% PARAMETERS
%
%   pos    - any position vector for packed array, e.g., G.faces.nodes.
%
%   offset - scalar or vector of rotation offsets. An offset of 1 give a
%            rotated index [n, 1, 2, ..., n-1] for each section.
%
% RETURNS
%
%   ix  - "rotated" vector of indices into value array, e.g., G.face.nodes
%
%
   pos    = int32(pos); 
   num    = diff(pos);
   offset = mod(int32(offset), num); % net offset
   ix     = zeros(max(pos)-1, 1);
   
   ix(mcolon(pos(1:end-1), pos(1:end-1)+num-offset-1)) = ...
      mcolon(pos(1:end-1)+offset, pos(2:end)-1);
   ix(mcolon(pos(1:end-1)+num-offset, pos(2:end)-1)) = ...
      mcolon(pos(1:end-1), pos(2:end)-1-num+offset);
end


