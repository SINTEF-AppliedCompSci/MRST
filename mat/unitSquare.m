function G = unitSquare(nx, ny)

[x, y] = meshgrid(0:1/(nx-1):1, 0:1/(ny-1):1);

x(:,2:nx-1) = x(:,2:nx-1) + random('Normal', 0, 1/(3*(nx-1)), nx, nx-2);
y(2:ny-1,:) = y(2:ny-1,:) + random('Normal', 0, 1/(3*(ny-1)), ny-2, ny);

P = [x(:), y(:)];
t = delaunayn(P);

G = triangleGrid(P,t);
G = pebi(G);

% X = G.nodes.coords;
% 
% Nn = G.nodes.num
% Ne = G.faces.num;
% 
% boundaryNodes = zeros(Nn,1);
% neighbors = G.faces.neighbors;
% for k = 1:Ne
%     if neighbors(k,1) == 0 || neighbors(k,2) == 0
%         nodeNum = G.faces.nodePos(k):G.faces.nodePos(k+1)-1;
%         boundaryNodes([G.faces.nodes(nodeNum)']) = 1;
%     end
% end
% 
% remNodes = (sum((X == 0) + (X == 1),2) > 0) - boundaryNodes;
% remNodes = find(remNodes);
% 
% G.nodes.coords(remNodes,:) = 4;

% G = removeNodes(G, remNodes);
% 
% function H = removeNodes(G, nodes)
%    map = false(G.nodes.num, 1);
%    map(nodes) = true;
% 
%    faceno = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .';
%    faces  = accumarray(faceno, map(G.faces.nodes), [G.faces.num,1],@(x) any(x));
% 
%    cellno = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
%    edges = double(reshape(G.faces.nodes,2,[])');
%    hf    = faces(G.cells.faces(:,1));


% 
%    [cellno(hf), map(edges(G.cells.faces(hf,1),:)), edges(G.cells.faces(hf,1),:)];
%    c        = repmat(cellno(hf), [1,2])';
%    e        = edges(G.cells.faces(hf,1),:)';
%    nodemask = map(e);
%    G = removeFaces(G, find(faces));
% 
%    tmp = [reshape(c(~nodemask), 2, [])',reshape(e(~nodemask), 2, [])'];
% 
%    neigh  = [tmp(:,1), zeros(size(tmp(:,1)))];
%    fnodes = reshape(tmp(:,3:4)', [], 1);
%    nnodes = repmat(2, [size(tmp, 1), 1]);
%    H = addFaces(G, fnodes, nnodes, neigh);
% 
%    H = removeCells(H, find(diff(H.cells.facePos)<3));
% end



end
