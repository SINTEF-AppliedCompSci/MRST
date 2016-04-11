function G = unitSquare(gridDim)

nx = gridDim(1); ny = gridDim(2);
dx = 1/nx; dy = 1/ny;

X = zeros(nx*ny,2);
yEdge = (dy/2:dy:(1-dy/2))';
xEdge = (dx/2:dx:(1-dx/2))';
yEdge = [repmat(yEdge,2,1); dy/2*ones(nx,1); (1-dy/2)*ones(nx,1)];
xEdge = [dx/2*ones(ny,1); (1-dx/2)*ones(ny,1); repmat(xEdge,2,1)];

x = [xEdge;rand(nx*ny-2*(nx+ny-2),1)*(1-2*dx) + dx];
y = [yEdge;rand(nx*ny-2*(nx+ny-2),1)*(1-2*dy) + dy];

X = [x,y];

t = delaunayn(X);
G = triangleGrid(X,t);
G = pebi(G);

end
