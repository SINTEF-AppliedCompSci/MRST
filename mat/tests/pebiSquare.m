function G = pebiSquare(gridDim, gridLim)
%   Generates a 2D pebi grid of gridDim(1) x gridDim(2) cells, covering the
%   rectangle [0, gridLim(1)] x [0, girdLim(2)].
%
%   SYNOPSIS:
%       G = unitSquar(gridDim, gridLim)
%
%   DESCRIPTION:
%       Generates 2D pebi grid with gridDim(1) x gridDim(2) cells covering
%       the rectangle [0, gridLim(1)] x [0, gridLim(2)] by setting equally
%       spaces point along the boundary, and random point inside the
%       domain.
%
%   REQUIRED PARAMETERS:
%       gridDim - vector of two elements specifying roughly number of
%       points ineach coordinate direction.
%       gridLim - Domain boundary.
%
%   RETURNS:
%       G   - MRST grid.  
%-----------------------------------------------------------------ØSK-2016-

nx = gridDim(1); ny = gridDim(2);
xMax = gridLim(1); yMax = gridLim(2);
dx = xMax/(nx-1); dy = yMax/(ny-1);

yEdge = (0:dy:yMax)';
xEdge = (0:dx:xMax)';
yEdge = [repmat(yEdge,2,1); zeros(nx,1); yMax*ones(nx,1)];
xEdge = [zeros(ny,1); xMax*ones(ny,1); repmat(xEdge,2,1)];

x = [xEdge;rand(nx*ny-2*(nx+ny-2),1)*(xMax-dx) + dx/2];
y = [yEdge;rand(nx*ny-2*(nx+ny-2),1)*(yMax-dy) + dy/2];

X = [x,y];
X = unique(X,'rows');

t = delaunayn(X);
G = triangleGrid(X,t);
G = pebi(G);

end
