%% Basic Grid Operations and Manipulations: Rectilinear/Curvilinear Grids
% In this tutorial, we describe how you can use MRST to construct various kinds of unstructured grids based on Delaunay triangulations and Voronoi diagrams. This example is discussed in more detail in 
% <http://www.sintef.no/projectweb/mrst/jolts/ Jolt 2>

%% Triangular grid in 2D
clf
[x,y] = meshgrid(0:10,0:8);
x(2:8,2:10) = x(2:8,2:10) + .25*randn(7,9);
y(2:8,2:10) = y(2:8,2:10) + .25*randn(7,9);

G = triangleGrid([x(:) y(:)]);
h=plotGrid(G,'FaceColor','none');

%% Voronoi diagram
V = pebi(G);
hold on, plot(x(:),y(:),'k.','MarkerSize',12); hold off
plotGrid(V,'FaceColor','none','LineWidth',2);
delete(h)


%% Standard data set from Matlab
load seamount;
plot(x(:),y(:),'k.','MarkerSize',12);

t = delaunay(x(:), y(:));
G = triangleGrid([x(:) y(:)], t);

plotGrid(G,'FaceColor','none');


%% Extrude to 3D
G = makeLayeredGrid(G,5);
clf
plotGrid(G,'FaceColor',[1 .8 .8]); view(5,45); axis tight off


%% Radial grid
P = [];
for r = exp(-3.5:0.25:0),
   [x,y,z] = cylinder(r,16); P = [P [x(1,:); y(1,:)]];          %#ok<AGROW>
end
[x,y] = meshgrid(-2.5:2.5,-2.5:2.5);
P = unique([P'; x(:) y(:); 0 0],'rows');

G = makeLayeredGrid(pebi(triangleGrid(P)), 3);
cla, plotGrid(G,'FaceColor',[1 .8 .8]); view(30,60), axis tight off

%% Cartesian grid with radial refinement
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
cla, plotGrid(G,'FaceColor','none'); axis([-.1 14.1 -.1 8.1]); axis off

%% Copyright notice

displayEndOfDemoMessage(mfilename)

% #COPYRIGHT_EXAMPLE#