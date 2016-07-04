%% Example 1: Cartesian grid with refinement in the middle
clf
G1 = cartGrid([ 4  4],[1 1]);
G2 = cartGrid([ 8  8],[1 1]);
G3 = cartGrid([12  4],[3 1]);

% Glue middle part
G1 = translateGrid(G1,[0 1]); plotGrid(G1,'FaceColor',[1 .9 .9]);
G2 = translateGrid(G2,[1 1]); plotGrid(G2,'FaceColor',[.9 1 .9]);
G = glue2DGrid(G1, G2);
G1 = translateGrid(G1,[2 0]); plotGrid(G1,'FaceColor',[1 .9 .9]);
G = glue2DGrid(G, G1);

% Glue top and bottom
G = glue2DGrid(G, G3);        plotGrid(G3,'FaceColor',[.9 .9 1]);
G3 = translateGrid(G3,[0 2]); plotGrid(G3,'FaceColor',[.9 .9 1]);
G = glue2DGrid(G3, G); %#ok<NASGU>
%print -depsc2 multiBlock-1a.eps;

%%
% Repeat the whole process without intermediate plotting
G1 = cartGrid([ 5  5],[1 1]);
G2 = cartGrid([20 20],[1 1]);
G3 = cartGrid([15  5],[3 1]);

G = glue2DGrid(G1, translateGrid(G2,[1 0]));
G = glue2DGrid(G,  translateGrid(G1,[2 0]));
G = glue2DGrid(G3, translateGrid(G, [0 1]));
G = glue2DGrid(G,  translateGrid(G3,[0 2]));
G = twister(G);
clf
plotGrid(G,'FaceColor','none');
%print -depsc2 multiBlock-1b.eps;

%% Example 2: Cartesian grid with triangular/PEBI refinement
% Construct unstructured grid
[N,M]=deal(10,15);
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

% Add cell tags
G2.cells.faces(:,2) = zeros(size(hf));
i = sum(N.*n,2)==-1; G2.cells.faces(i,2) = 1;
i = sum(N.*n,2)== 1; G2.cells.faces(i,2) = 2;
n = n(:,[2 1]);
i = sum(N.*n,2)==-1; G2.cells.faces(i,2) = 3;
i = sum(N.*n,2)== 1; G2.cells.faces(i,2) = 4;

% Glue grids together
G = glue2DGrid(G1, translateGrid(G2,[1 0]));
G = glue2DGrid(G,  translateGrid(G1,[2 0]));
G = glue2DGrid(G3, translateGrid(G, [0 1]));
G = glue2DGrid(G,  translateGrid(G3,[0 2]));
G = twister(G);
clf
plotGrid(G,'FaceColor','none');
%print -depsc2 multiBlock-2a.eps;

%% 
% Repeat with Voronoi grid instead
[N,M]=deal(8,12);
xv = linspace(0,1,N+1);
yv = linspace(0,1,M+1);
[x,y] = ndgrid(xv, yv);
x(2:N,2:M) = x(2:N,2:M) + 0.3*randn(N-1,M-1)*max(diff(xv));
y(2:N,2:M) = y(2:N,2:M) + 0.3*randn(N-1,M-1)*max(diff(yv));
G2 = computeGeometry(pebi(triangleGrid([x(:) y(:)])));

% Set orientation of faces
hf = G2.cells.faces(:,1);
hf2cn = gridCellNo(G2);
sgn = 2*(hf2cn == G2.faces.neighbors(hf, 1)) - 1;
N   = bsxfun(@times, sgn, G2.faces.normals(hf,:));
N   = bsxfun(@rdivide, N, G2.faces.areas(hf,:));
n   = zeros(numel(hf),2); n(:,1)=1;

G2.cells.faces(:,2) = zeros(size(hf));
i = sum(N.*n,2)==-1; G2.cells.faces(i,2) = 1;
i = sum(N.*n,2)== 1; G2.cells.faces(i,2) = 2;
n = n(:,[2 1]);
i = sum(N.*n,2)==-1; G2.cells.faces(i,2) = 3;
i = sum(N.*n,2)== 1; G2.cells.faces(i,2) = 4;

% Refine G3 in the x-direction and rescale in y-direction
G3 = cartGrid([24 5], [3 .25]);
G = glue2DGrid(G1, translateGrid(G2,[1 0]));
G = glue2DGrid(G,  translateGrid(G1,[2 0]));
G = glue2DGrid(G3, translateGrid(G, [0 .25]));
G = glue2DGrid(G,  translateGrid(G3,[0 1.25]));
G = twister(G);
clf
plotGrid(G,'FaceColor','none');
%print -depsc2 multiBlock-2b.eps;

%% Example 3: Extruded grid rotated
G = glue2DGrid(G1, translateGrid(G2,[0 1]));
G = glue2DGrid(G,  translateGrid(G1,[0 2]));
G = makeLayeredGrid(G, 5);
G.nodes.coords = G.nodes.coords(:,[3 1 2]);
clf, plotGrid(G,'FaceColor',[1 1 1]); view(115,20); axis off
%print -depsc2 multiBlock-3.eps;
