function G = sortEdges(G)

%{ 
Copyright 2009-2014 SINTEF ICT, Applied Mathematics
%} 
assert(G.griddim==2);%??? is this correct
G=sortCellFaces(G);
end
%% Helper function to sortCellFaces
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

%% For each (2D) cell, sort edges counter-clockwise to ensure that the cell
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

   for i=1:G.cells.num%size(facePos, 1)-1,
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
      cellFaces(I,1) = cf;
      if(size(cellFaces,2)==2)
        cellFaces(I,2) = cellFaces(I(m),2);
      end
   end

   G.cells.faces = cellFaces;
   G.faces.nodes = reshape(Edges', [], 1);
end
