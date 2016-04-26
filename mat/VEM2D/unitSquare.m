function G = unitSquare(gridDim, gridLim)

nx = gridDim(1); ny = gridDim(2);
xMax = gridLim(1); yMax = gridLim(2);
dx = xMax/(nx-1); dy = yMax/(ny-1);

% yEdge = (dy/2:dy:(yMax-dy/2))';
yEdge = (0:dy:yMax)';
% xEdge = (dx/2:dx:(xMax-dx/2))';
xEdge = (0:dx:xMax)';
% yEdge = [repmat(yEdge,2,1); dy/2*ones(nx,1); (yMax-dy/2)*ones(nx,1)];
yEdge = [repmat(yEdge,2,1); zeros(nx,1); yMax*ones(nx,1)];
% xEdge = [dx/2*ones(ny,1); (xMax-dx/2)*ones(ny,1); repmat(xEdge,2,1)];
xEdge = [zeros(ny,1); xMax*ones(ny,1); repmat(xEdge,2,1)];

x = [xEdge;rand(nx*ny-2*(nx+ny-2),1)*(xMax-dx) + dx/2];
y = [yEdge;rand(nx*ny-2*(nx+ny-2),1)*(yMax-dy) + dy/2];

X = [x,y];
X = unique(X,'rows');

t = delaunayn(X);
G = triangleGrid(X,t);
G = pebi(G);

end
