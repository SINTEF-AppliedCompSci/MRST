%% Grids from standard data sets in MATLAB
load trimesh2d;
G = triangleGrid([xfe yfe], trife); figure, plotGrid(G);

load tetmesh;
G = tetrahedralGrid(X,tet); figure; plotGrid(G); view(3);

%% Extend triangleGrid to triangulated surfaces
load trimesh3d;
G = myTriangleGrid([x y z], tri); plotGrid(G); view(3);

%% Local grid refinement (almost structured)
[x y] = meshgrid(-.25:.5:6.25);
pt = [x(:) y(:)];
i=all((pt<5) & (pt>1),2);
pt = pt(~i,:);
[X Y] = meshgrid(1.125:.25:4.875);
pt = unique([pt; X(:) Y(:)],'rows');
i=all((pt<4) & (pt>2),2);
pt = pt(~i,:);
[X Y] = meshgrid(2.0625:.125:3.9375);
pt = unique([pt; X(:) Y(:)],'rows');
G = pebi(triangleGrid(pt));
clf
G = removeCells(G, unique(sum(G.faces.neighbors(any(G.faces.neighbors==0,2),:),2)));
plotGrid(G, 'FaceColor','none','LineWidth',2); axis tight

%% Padded, hexagonal grid adapted to two faults
[x y] = meshgrid(-0.25:.5:8.75);
pt = [x(:) y(:)];
i=all((pt(:,1)<8) & (pt(:,2)<7.5) & (pt(:,1)>1) & (pt(:,2)>1),2);
pt = pt(~i,:);
[x,y] = meshgrid((.75:.5:4.5)*2*cos(pi/6),(1:0.5:7)*2*sin(pi/6));
x = [x(:); x(:)+.5*cos(pi/6)];
y = [y(:); y(:)+.5*sin(pi/6)];
pt = [pt; x(:) y(:)];
x1 = cos(pi/6)*(2.5:.25:4.5).';
y1 = 2+sin(pi/6)*(0:.25:2).';
d1 = distancePt2Line(pt',repmat([x1(end) y1(end)]',1,size(pt,1)), ...
                        repmat([x1(1) y1(1)]'+eps,1,size(pt,1))).';
x2 = .5*cos(pi/2)+cos(pi/6)*(2.5:.25:7.5).';
y2 = 5+.5*sin(pi/2)-sin(pi/6)*(0:.25:5).';
d2 = distancePt2Line(pt',repmat([x2(end) y2(end)]',1,size(pt,1)), ...
                        repmat([x2(1) y2(1)]'+eps,1,size(pt,1))).';
pt = vertcat(pt((d1>.45) & (d2>.45),:), ...
             [x1+.125*cos( 2*pi/3) y1+.125*sin( 2*pi/3)], ...
             [x1+.125*cos(  -pi/3) y1+.125*sin(  -pi/3)], ...
             [x1+.375*cos( 2*pi/3) y1+.375*sin( 2*pi/3)], ...
             [x1+.375*cos(  -pi/3) y1+.375*sin(  -pi/3)], ...
             [x2+.125*cos(   pi/3) y2+.125*sin(   pi/3)], ...
             [x2+.125*cos(-2*pi/3) y2+.125*sin(-2*pi/3)], ...
             [x2+.375*cos(   pi/3) y2+.375*sin(   pi/3)], ...
             [x2+.375*cos(-2*pi/3) y2+.375*sin(-2*pi/3)]);
pt = unique(pt,'rows');
G = pebi(triangleGrid(pt));
clf
G = removeCells(G, unique(sum(G.faces.neighbors(any(G.faces.neighbors==0,2),:),2)));
plotGrid(G, 'FaceColor','none','LineWidth',1); axis tight
%hold on; plot(pt(:,1),pt(:,2),'ro'); hold off
hold on; plot(x1,y1,'-k',x2,y2,'-k','LineWidth',3); hold off;