function G = unitSquare(nx, ny)

dx = 1/(nx-1); dy = 1/(ny-1);
[x, y] = meshgrid(-dx/2:dx:1+dx/2, -dy/2:dy:1+dy/2);

x(3:ny-1,3:nx-1) = x(3:ny-1,3:nx-1) + random('Normal', 0, dx/4, ny-3, nx-3);
y(3:ny-1,3:nx-1) = y(3:ny-1,3:nx-1) + random('Normal', 0, dy/4, ny-3, nx-3);

P = [x(:), y(:)];
t = delaunayn(P);

G = triangleGrid(P,t);
G = pebi(G);

neighbors = G.faces.neighbors;
tmp = (neighbors(:,1) == 0) + (neighbors(:,2) == 0);
remCells = unique(sum(neighbors(find(tmp),:),2 ));
G = removeCells(G,remCells);

%   distmesh

% X = G.nodes.coords;
% G.nodes.coords(find(X(:,1) < 0),1) = 0;
% G.nodes.coords(find(X(:,1) > 1),1) = 1;
% G.nodes.coords(find(X(:,2) < 0),2) = 0;
% G.nodes.coords(find(X(:,2) > 1),2) = 1;

end
