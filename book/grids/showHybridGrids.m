%% Hybrid grids

%% Example 1: Cartesian grid with radial refinement
Pw = [];
for r = exp(-3.5:.2:0),
    [x,y,z] = cylinder(r,28); Pw = [Pw [x(1,:); y(1,:)]];
end
Pw = [Pw [0; 0]];
Pw1 = bsxfun(@plus, Pw, [2; 2]);
Pw2 = bsxfun(@plus, Pw, [12; 6]);
[x,y] = meshgrid(0:.5:14, 0:.5:8);
P = unique([Pw1'; Pw2'; x(:) y(:)], 'rows');
G = pebi(triangleGrid(P));

%% Example 2: Glue multiple grids together
% Construct unstructured grid
[N,M]=deal(20,10);
xv = linspace(0,1,N+1);
yv = linspace(0,1,M+1);
[x,y] = ndgrid(xv, yv);
x(2:N,2:M) = x(2:N,2:M) + 0.3*randn(N-1,M-1)*max(diff(xv));
y(2:N,2:M) = y(2:N,2:M) + 0.3*randn(N-1,M-1)*max(diff(yv));
G1 = computeGeometry(triangleGrid([x(:) y(:)]));

% Set orientation of faces
hf = G1.cells.faces(:,1);
hf2cn = gridCellNo(G1);
sgn = 2*(hf2cn == G1.faces.neighbors(hf, 1)) - 1;
N   = bsxfun(@times, sgn, G1.faces.normals(hf,:));
N   = bsxfun(@rdivide, N, G1.faces.areas(hf,:));
n   = zeros(numel(hf),2); n(:,1)=1;

G1.cells.faces(:,2) = zeros(size(hf));
i = sum(N.*n,2)==-1; G1.cells.faces(i,2) = 1; sum(i)
i = sum(N.*n,2)== 1; G1.cells.faces(i,2) = 2; sum(i)
n = n(:,[2 1]);
i = sum(N.*n,2)==-1; G1.cells.faces(i,2) = 3; sum(i)
i = sum(N.*n,2)== 1; G1.cells.faces(i,2) = 4; sum(i)

% Second grid
[N,M]=deal(15,8);
xv = linspace(0,1,N+1);
yv = linspace(0,1,M+1);
[x,y] = ndgrid(xv, yv);
x(2:N,2:M) = x(2:N,2:M) + 0.3*randn(N-1,M-1)*max(diff(xv));
y(2:N,2:M) = y(2:N,2:M) + 0.3*randn(N-1,M-1)*max(diff(yv));
G2 = computeGeometry(triangleGrid([x(:) y(:)]));

% Set orientation of faces
hf = G2.cells.faces(:,1);
hf2cn = gridCellNo(G2);
sgn = 2*(hf2cn == G2.faces.neighbors(hf, 1)) - 1;
N   = bsxfun(@times, sgn, G2.faces.normals(hf,:));
N   = bsxfun(@rdivide, N, G2.faces.areas(hf,:));
n   = zeros(numel(hf),2); n(:,1)=1;

G2.cells.faces(:,2) = zeros(size(hf));
i = sum(N.*n,2)==-1; G2.cells.faces(i,2) = 1; sum(i)
i = sum(N.*n,2)== 1; G2.cells.faces(i,2) = 2; sum(i)
n = n(:,[2 1]);
i = sum(N.*n,2)==-1; G2.cells.faces(i,2) = 3; sum(i)
i = sum(N.*n,2)== 1; G2.cells.faces(i,2) = 4; sum(i)


% Construct Cartesian grid
G3 = computeGeometry(cartGrid([25 8],[1 1]));
G4 = computeGeometry(cartGrid([14 7],[1 1]));



%%
v = [0 1; 0 -1; -1 0; 1 0];
clf
for i=1:4
    g1 = G1;
    g2 = G2;
    g1.nodes.coords = bsxfun(@plus,g1.nodes.coords,v(i,:));
    g1 = computeGeometry(g1);
    subplot(1,2,1); cla
    plotGrid(g1,'FaceColor',[.7 1 .7]);
    plotGrid(g2,'FaceColor',[.7 .7 1]); axis equal off
    G = glue2DGrid(g2,g1);
    subplot(1,2,2); cla;
    plotGrid(G); axis equal off
    drawnow; pause;
end
