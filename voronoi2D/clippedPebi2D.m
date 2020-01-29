function G = clippedPebi2D(p, bnd)
% Construct a 2D clipped Pebi grid.
%
% SYNOPSIS:
%   G = clippedPebi2D(p, bnd)
%
% PARAMETERS:
%   p      - A n X 2 array of coordinates. Each coordinate coresponds to 
%            one Voronoi site.
%   bnd    - A k X 2 array of coordinates. Each coordinate coresponds to a
%            vertex in the polygon boundary. The coordinates must be
%            ordered clockwise or counter clockwise. 
%
% RETURNS:
%   G      - Valid MRST grid definition.
%
% EXAMPLE:
%   p = rand(30,2);
%   bnd = [0,0;0,1;1,1;1,0];
%   G = clippedPebi2D(p,bnd);
%   plotGrid(G);
%   axis equal tight
%
% SEE ALSO
%   compositePebiGrid2D, pebi, pebiGrid2D, clippedPebi3D

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2015-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

dt = delaunayTriangulation(p);
E = dt.edges();

m = size(bnd,1);



rem = true(size(E,1),1);
C = cell(size(p,1),1);
V = [];
for s = 1:size(p,1)
    NC = [E(:,2)==s, E(:,1)==s];
    bisect = find(any(NC,2));
  
    n = bsxfun(@minus, p(E(NC),:), p(s,:));
    n = bsxfun(@rdivide, n,sqrt(sum(n.^2,2)));
    x0 =bsxfun(@plus, p(E(NC),:), p(s,:))/2;
    
    symT = repmat({-ones(1,3)},size(bnd,1),1);

    [newVertex, ~] = clipPolygon(bnd, n, x0, symT, bisect,'noSym',true);

    newVertex = round(newVertex*1e6)/1e6;
    newVertex = unique(newVertex,'rows','stable');
    %keep = inpolygon(newVertex(:,1), newVertex(:,2), bnd(:,1), bnd(:,2));
    %newVertex = newVertex(keep,:); % The clipPolygon routine might return
    if isempty(newVertex)          % vertexes outside the domain. The part
      rem(any(NC,2)) = false;      % of the cell outside the domain will have
      continue                     % zero volume. 
    end
    C{s} = [C{s}, size(V,1)+1:size(V,1)+size(newVertex,1)];
    V = [V;newVertex];
    

end
V = round(V*10^6)/10^6;
[V,~,IC] = unique(V,'rows');

G.cells.num = numel(C);
G.cells.facePos = cumsum([1; cellfun(@numel, C)]);
C2 = cellfun(@(c) [c,c(1)],C,'un', false);
faces = horzcat(C2{:})';
faces = [faces(1:end-1), faces(2:end)];
addTo = cumsum([0;ones(G.cells.num-2,1)]);

faces(G.cells.facePos(2:end-1) + addTo,:) = [];
faces = sort(IC(faces),2);

[nodes,~,faces] = unique(faces, 'rows');

G.cells.faces = faces;

G.nodes.coords = V;
G.nodes.num = size(V,1);

G.faces.num = max(faces);
G.faces.nodes = reshape(nodes',[],1);
G.faces.nodePos = (1:2:2*G.faces.num+1)';


cellNo            = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';
G.faces.neighbors = zeros(G.faces.num,2);
for i = 1:G.faces.num
    neigh = G.cells.faces==i;
    if sum(neigh)==2
        G.faces.neighbors(i,[1,2]) = cellNo(neigh);
    else
          G.faces.neighbors(i,:) = [cellNo(neigh),0];
    end
end

G.griddim = 2;
G.type = {mfilename};
G = sortEdges(G);
end
