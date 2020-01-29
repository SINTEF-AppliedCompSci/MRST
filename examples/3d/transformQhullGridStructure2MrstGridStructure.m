%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%% Create voronoi Grid:
% We first create a Voronoi grid using qhull. 
n   = 5;
[X,Y,Z] = ndgrid(linspace(0,1,n));
pts = [X(:),Y(:),Z(:)];
pts(1:2:end,1) = pts(1:2:end,1) + 0.1;
[V,C] = voronoin(pts);

%% Remove infinity cells
% Our MRST grid structure does not support cells that extend to infinity.
% We therefore remove these before we convert the grid structure.
rem = cellfun(@(c) any(isinf(V(c,1))), C);
Cn   = C(~rem);
Vn   = V(2:end,:);
Cn   = cellfun(@(c) c-1, Cn,'un',false);
G1   = voronoi2mrstGrid3D(Vn, Cn);
% Or we can let voronoi2mrstGrid3D remove these. Note: this will result in a
% warning.
G2 = voronoi2mrstGrid3D(V,C);

%% Plot grid
clf;
plotGrid(G1)
axis equal tight
view(40,30)
light()
