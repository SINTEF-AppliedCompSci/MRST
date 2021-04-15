function H = pebi(G, varargin)
%Compute dual grid of triangular grid G.
%
% SYNOPSIS:
%   H = pebi(G)
%
% PARAMETERS:
%   G       - Triangular grid structure as described by grid_structure.
%
% RETURNS:
%   H       - Dual of a triangular grid with edges that are perpendicular
%             bisectors of the edges in the triangular grid.  This grid is
%             also  known as a Voronoi diagram.
%
%
% NOTE:
%   If triangulation is not Delaunay, some circumcenters may fall outside
%   its assiciated triangle.  This may generate warped grids or grids that
%   do not preserve the original outer boundary of the triangulation.
%
% SEE ALSO:
%   `triangleGrid`

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

   p         = circumcenter(G);
   [G, p, a] = addAuxillaryTriangles(G, p);
   
   % Add midpoints of boundary edges in G
   b         = any(G.faces.neighbors==0, 2);
   fe        = reshape(G.faces.nodes, 2, []) .';

   % Midpoints of all G- boundary edges
   newcoords = 0.5*(G.nodes.coords(fe(b,1),:) + G.nodes.coords(fe(b,2),:));
   newnodes  = G.cells.num + (1:size(newcoords, 1))';
   newnodes_art = a(fe(b,1)) | a(fe(b,2));

   % Boundary nodes in G
   gcoords   = G.nodes.coords(unique(fe(b,:)),:);
   gcoords_a = a(unique(fe(b,:)));
   map       = cumsum(accumarray(reshape(fe(b,:),[], 1), 1)>0);
   num       = size(p,1)+size(newcoords, 1);

   % Edges in H:
   % internal edges
   faceNodes = G.faces.neighbors;

   % G.faces.neigbors == edges in H.  Entries with 0 (on boundary) has no
   % second node.  Add midpoint of G-boundary nodes.
   faceNodes(b, :) = [sum(faceNodes(b,:), 2), newnodes];

   z = zeros(sum(b), 1);
   neighbors = [reshape(G.faces.nodes, 2, [])';...G-nodes == H-cells
                z,       fe(b,1); ...
                fe(b,2), z      ];
   faceNodes = [faceNodes ;...
                num+map(fe(b,1)), newnodes;...
                num+map(fe(b,2)), newnodes];
   faceNodes = reshape(faceNodes', [], 1);

   % Wrap up
   numFaces  = size(neighbors, 1);
   numCells  = max(neighbors(:));
   numNodes  = repmat(2, [numFaces, 1]);
   H         = emptyGrid(numCells, [p; newcoords;gcoords]);
   H         = addFaces(H,faceNodes, numNodes, neighbors);

   H = sortCellFaces(H);

   % Fix boundary
   a = [false(size(p, 1), 1);newnodes_art; gcoords_a];
   % Remove nodes 'a' in H
   H = removeNodes(H, find(a));


   % Make cells oriented
   assert(H.nodes.num==numel(unique(H.faces.nodes)));
   assert(H.nodes.num==max(H.faces.nodes));
   H = sortCellFaces(H);
   H = mergeBoundaryFaces(H);
   assert(H.nodes.num==numel(unique(H.faces.nodes)));
   assert(H.nodes.num==max(H.faces.nodes));
   % H = sortCellFaces(H);
   H.type    = { mfilename };
   H.griddim = 2;

   H = removeDuplicateNodes(H);
   
   H=sortCellFaces(H);
   assert(H.nodes.num==numel(unique(H.faces.nodes)));
   assert(H.nodes.num==max(H.faces.nodes));
end

function H = removeDuplicateNodes(H)
% Fix node uniqueness.  This is important since we use node uniqueness to
% determine connectivity of faces.

   assert(H.nodes.num == size(H.nodes.coords, 1));
   [nodes, map, map] = unique(H.nodes.coords, 'rows');                     %#ok<ASGLU>

   if size(nodes, 1) < H.nodes.num,
       % map facenodes to new node numbers
       assert(all(all(nodes(map,:) == H.nodes.coords)));
       H.nodes.coords     = nodes;
       H.nodes.num        = size(nodes, 1);
       H.faces.nodes(:,1) = map(H.faces.nodes(:,1));

       % remove collapsed faces
       f = find(H.faces.nodes(1:2:end)==H.faces.nodes(2:2:end));
       H = removeFaces(H, f);

       % remove cells with fewer than 3 faces
       c = diff(H.cells.facePos) < 3;
       H = removeCells(H, c);
   end
end

function [G, p, a] = addAuxillaryTriangles(G, p)
   G = sortEdges(G);
   t             = grid2tri(G);
   [inside, col] = pointInside(t, G.nodes.coords, p);
   % Compute half-face from inside and col. (make sure col is suitable)
   % If circumcenter falls outside domain, make sure boundary remains
   % intact by adding auxillary triangles.
   hf            = 3*find(~inside)-3+mod(col,3)+1;
   isbd          = any(G.faces.neighbors(G.cells.faces(hf,1), :)==0, 2);

   % Fix only points that are on boundary of domain to ensure boundary
   % remains intact for dual grid.
   a = false(G.nodes.num, 1);
   if any(isbd),
      f = find(~inside);
      f = f(isbd);
      col = col(isbd);
      ix1 =t(sub2ind(size(t), f, col));
      ix2 =t(sub2ind(size(t), f, mod(col,3)+1));
      ix3 =t(sub2ind(size(t), f, mod(col+1,3)+1));

      x   = G.nodes.coords(:,1);
      y   = G.nodes.coords(:,2);
      u   = [x(ix2)-x(ix1), y(ix2)-y(ix1)];
      v   = [x(ix3)-x(ix1), y(ix3)-y(ix1)];
      w   = v - bsxfun(@times, sum(u.*v, 2)./sum(u.^2,2), u);

      p = [x(ix3),y(ix3)]-2*w;

      num = G.nodes.num+(1:sum(isbd));
      x   = [G.nodes.coords(:,1);p(:, 1)];
      y   = [G.nodes.coords(:,2);p(:, 2)];
      T   = t(f,:);T = [T, T];
      % >For loop
      i=(1:numel(col))';
      T(sub2ind(size(T),i,col(i)))=num(i);
      T(sub2ind(size(T),i,3+mod(col(i),3)+1))=num(i);

      T = reshape(T', 3, [])';
      i=false(size(t,1),1);i(f)=true;
      t = [t(~i, :); T];
      clear G
      G = triangleGrid([x,y], t);      % New trigrid with auxillary triangles on bd.
      a(num) = true;
      p = circumcenter(G);  % Circumcenters now fall on boundary.
   end
end


% For cells with two boundary faces, check if faces are parallel... if so,
% replace them with a single face.
function G = mergeBoundaryFaces(G)
   cellno = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   sgn    = (cellno == G.faces.neighbors(G.cells.faces(:,1), 1));

   b  = any(G.faces.neighbors == 0, 2);
   i  = find(b(G.cells.faces(:,1)));
   f  = G.cells.faces(i,1);
   c  = sum(G.faces.neighbors(f,:),2);
   e  = reshape(G.faces.nodes, 2, [])';
   be = e(f,:); be(~sgn(i), :) = be(~sgn(i), [2,1]);

   % Pick out faces (cells and edges) with two boundary faces.
   ind = c(1:end-1) == c(2:end);
   ind = any([[ind; true], [true;ind]], 2);
   E = be(ind, :);
   f = f(ind);
   c = c(ind);

   % Remove those pairs whose faces are parallel...
   isparallel = compareNormals(E(1:2:end-1,:), E(2:2:end,:), ...
                               G.nodes.coords);
   f = f(rldecode(isparallel, 2));
   H = removeFaces(G, f);

   % ...and replace with one face.
   fnodes = reshape(E(rldecode(isparallel, 2),:)', 4, [])';
   fnodes = reshape(fnodes(:,[1,4])',[],1);
   nnodes = repmat(2, [numel(fnodes)/2, 1]);
   c = c(1:2:end);
   neigh  = zeros(sum(isparallel), 2);
   neigh(:,1) = c(isparallel);

   G = addFaces(H, fnodes, nnodes, neigh);
   G=remapNodes(G);
end

% Remove nodes N.  Remove faces and cells associated with N.
function H = removeNodes(G, nodes)
   map = false(G.nodes.num, 1);
   map(nodes) = true;

   faceno = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .';
   faces  = accumarray(faceno, map(G.faces.nodes), [G.faces.num,1],@(x) any(x));

   cellno = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   edges = double(reshape(G.faces.nodes,2,[])');
   hf    = faces(G.cells.faces(:,1));

   %[cellno(hf), map(edges(G.cells.faces(hf,1),:)), edges(G.cells.faces(hf,1),:)];
   c        = repmat(cellno(hf), [1,2])';
   e        = edges(G.cells.faces(hf,1),:)';
   nodemask = map(e);
   G = removeFaces(G, find(faces));

   tmp = [reshape(c(~nodemask), 2, [])',reshape(e(~nodemask), 2, [])'];

   neigh  = [tmp(:,1), zeros(size(tmp(:,1)))];
   fnodes = reshape(tmp(:,3:4)', [], 1);
   nnodes = repmat(2, [size(tmp, 1), 1]);
   H = addFaces(G, fnodes, nnodes, neigh);
   H = removeCells(H, find(diff(H.cells.facePos)<3));
   % renumber nodes all nodes has to be on a face
   %{
   [newnodes,ia,ic]=unique(G.faces.nodes);
   %
   if(numel(newnodes)< G.nodes.num)
    G.faces.nodes=G.faces.node(ic);
    G.nodes.coords=G.nodes.coords(newnodes,:);
    G.node.num=numel(newnodes);
   end
   %}
   H=remapNodes(H);   
end

function G=remapNodes(G)

   [newnodes,ia,ic]=unique(G.faces.nodes);                                 %#ok<ASGLU>
   %
   if(numel(newnodes)< G.nodes.num)
    G.faces.nodes=ic;
    G.nodes.coords=G.nodes.coords(newnodes,:);
    G.nodes.num=numel(newnodes);
   end
   assert(G.nodes.num==numel(unique(G.faces.nodes)));
   assert(G.nodes.num==max(G.faces.nodes));
end

function t = grid2tri(G)
% Extract (n x 3) array of triangle node numbers from grid
   assert(all(diff(G.faces.nodePos) == 2)); % 2D-grid
   assert(all(diff(G.cells.facePos) == 3)); % 3 edges per cell.

   cellno   = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   reverse  = (cellno == G.faces.neighbors(G.cells.faces(:,1), 2));
   edges     = reshape(G.faces.nodes, 2, []) .'; % Each row is a face edge.
   e         = edges(G.cells.faces(:,1), [2,1]);
   e(reverse,:) = e(reverse,[2,1]);
   t = reshape(e(:,1), 3, [])';
end

function isparallel = compareNormals(nn1,nn2, coord)
% Check if normals of edges nn1 and nn2 are parallel.  nn1 and nn2 are
% (n x 2) arrays of node numbers.
   x  = coord(:,1);
   y  = coord(:,2);
   n  = [diff(y(nn1), 1, 2), -diff(x(nn1), 1, 2)];
   e  = [diff(x(nn2), 1, 2),  diff(y(nn2), 1, 2)];
   ip = abs(sum(n.*e, 2));
   isparallel = (ip < 10*sqrt(eps));
end

function p = circumcenter(G)
% Find cell circumcenter, i.e., center of circle passing through each node
% in (triangular) cell.
   t     = reshape(G.cells.faces(:,1), 3, [])';
   edges = reshape(G.faces.nodes, 2, [])';
   E      = edges(t',:);
   cellno = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   reverse = cellno == G.faces.neighbors(G.cells.faces(:,1), 2);
   E(reverse, :) = E(reverse, [2,1]);
   tmp = reshape(E', 6, []);%tmp   = sort(reshape(E', 6, []));
   trinodes = tmp(1:2:5, :)';

   B = G.nodes.coords(trinodes(:,2), :) - G.nodes.coords(trinodes(:,1),:);
   C = G.nodes.coords(trinodes(:,3), :) - G.nodes.coords(trinodes(:,1),:);
   D = 2*(B(:,1).*C(:,2) - B(:,2).*C(:,1));
   p(:,1) = (C(:,2).*sum(B.^2, 2) - B(:,2).*sum(C.^2, 2))./D + G.nodes.coords(trinodes(:,1),1);
   p(:,2) = (B(:,1).*sum(C.^2, 2) - C(:,1).*sum(B.^2, 2))./D + G.nodes.coords(trinodes(:,1),2);

   %{
   % Check if circumcenter is inside triangle
   inside = pointInside(trinodes, G.nodes.coords, p);
   if ~all(inside),
      warning('Some circumcenters are outside corresponding triangle.');
   end
   %}
end

function [t, n] = pointInside(t, p, points)
% For each triangle t(k,:), check if point p(k,:) is inside t(k,:).
% p : (n x 2) are node coordinates for the triangle list t : (m x 3).
   x = p(:,1); y = p(:,2);
   vx = reshape(x(t)',3,[])'-repmat(points(:,1), [1,3]);
   vy = reshape(y(t)',3,[])'-repmat(points(:,2), [1,3]);

   wx = diff(reshape(x([t, t(:,1)])',4,[]), 1, 1)';
   wy = diff(reshape(y([t, t(:,1)])',4,[]), 1, 1)';

   a = (vx.*wy - vy.*wx);
   t = all(a > -sqrt(eps), 2);
   assert(all(sum(sign(a(~t,:)') == -1, 1) == 1), ...
      'Internal error in pebi.m');
   [n,dummy] = find(sign(a(~t,:)') == -1); %#ok

end


%{
function p = findCellPoints(G)
% Find average point in cell.
   d = size(G.nodes.coords, 2);

   e1 = G.faces.nodes(G.faces.nodePos(1:end-1));
   e2 = G.faces.nodes(G.faces.nodePos(1:end-1)+1);

   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   nodes = zeros(G.cells.num, d);
   for i=1:d
      nodes(:,i) = accumarray(cellNo, G.nodes.coords(e1(G.cells.faces(:,1),1), i))...
                 + accumarray(cellNo, G.nodes.coords(e2(G.cells.faces(:,1),1), i));
   end

   denominator = accumarray(cellNo, 1);
   p           = bsxfun(@rdivide, nodes, 2*denominator);
end
%}


% Helper function to sortCellFaces
function [edges, m, s] = swap(k, edges, m, s)
% Do one edge swap.  With N=size(edges,1), N-1 edge swaps will sort edges.

   % Find next edge
   i = k+find(any(edges(k+1:end,:)==k, 2));

   if any(i),
      i = i(1);

      % Swap rows i and k
      a=edges(k,:); edges(k,:)=edges(i,:); edges(i,:)=a;
      a=m(k);       m(k)=m(i);             m(i)=a;
      a=s(k);       s(k)=s(i);             s(i)=a;
   end

   % Flip direction of edge and corresponding sign.
   if edges(k,1)~=k,
      edges(k,:) = edges(k,[2,1]);
      s(k)    = -1;
   end

   % For algorithmic convenience, swap entries in node mapping (if we are
   % not at last edge)
   if edges(k,2)>k+1,
      u = edges == k+1;
      v = edges == edges(k,2);
      edges(u)  = edges(k,2);
      edges(v)  = k+1;
   end
end

%  For each (2D) cell, sort edges counter-clockwise to ensure that the cell
%  volumes are positive and that the edge normals and edge signs are
%  consistent with the cell-face neightbor list.
function G = sortCellFaces(G)
% Permute G.cells.faces and G.faces.nodes such that the faces are oriented
% counter clockwise when the face sign is used to switch edge direction.
%
% In a 2D grid, each face has two nodes.  The faces of each cell that is
% stored in G.cells.faces, is assumed to be sorted in the sense that the
% last node of one face is the first node of the next, e.g.,
%

   facePos   = G.cells.facePos;
   cellFaces = G.cells.faces;
   Edges     = reshape(G.faces.nodes, 2, []) .'; % Each row is a face edge.
   x         = G.nodes.coords(:,1);
   y         = G.nodes.coords(:,2);

   for i=1:size(facePos, 1)-1,
      I     = facePos(i):facePos(i+1)-1;
      cf    = cellFaces(I,1);

      % Find start and end node of each edge (face, that is) in a cell,
      % arranged as rows in an (numfaces x 2) array. Switch direction of
      % each edge (face) that has negative sign (in the sense given by
      % G.faces.neighbors)
      reverse           = (i == G.faces.neighbors(cf, 2));
      edges             = Edges(cf,:);
      edges(reverse, :) = edges(reverse, [2,1]);

      % Map node numbers to local node numbers 1:(num of different nodes)
      map        = accumarray(edges(:), 1);
      map(map>0) = (1:sum(map>0))';
      edges      = map(edges);

      % Sort edges e.
      m = 1:numel(cf);         % edge permutation
      s = ones(numel(cf),1);   % edge sign
      for k=1:numel(cf),
         [edges,m,s] = swap(k, edges, m, s);
      end
      cf = cf(m);

      % Sorted edges with global node numbers.
      assert(all((i==G.faces.neighbors(cf,2)) == reverse(m)))

      edges                = Edges(cf,:);
      edges(reverse(m), :) = edges(reverse(m), [2,1]);
      edges(s==-1,:)       = edges(s==-1, [2,1]);

      % Ensure that area is positive
      area    = 0.5* sum(sum(x(edges), 2).*diff(y(edges), 1, 2));
      if area<0,
         cf =  flipud(cf);
         s  =  flipud(-s);
      end

      % All edges in e should be oriented from columns 1 to 2, and
      % cell e should be oriented when parsed from first to last edge.
      assert( all(  edges(:, 2)==edges([2:numel(cf),1], 1)  ) )

      % Find edges that should be explicitly reversed, i.e., edges where
      % the sign s == -1:
      ind          = s==-1;
      Edges(cf(ind),:) = Edges(cf(ind),[2,1]);

      % Store result
      cellFaces(I) = cf;
   end

   G.cells.faces = cellFaces;
   G.faces.nodes = reshape(Edges', [], 1);
end
