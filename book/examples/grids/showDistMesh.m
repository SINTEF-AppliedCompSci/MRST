%% Square domain with circular cutout, uniform mesh
mrstModule add distmesh;
fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
[p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
G = triangleGrid(p, t);
clf, plotGrid(G); axis equal tight off;

%% Same domain, but graded grid
fh=@(p) 0.05+0.3*dcircle(p,0,0,0.5);
[p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
G = triangleGrid(p, t);
clf, plotGrid(G); axis equal tight off;

%% Polyhedral domain, graded PEBI grid
pv=[-1 -1; 0 -.5; 1 -1; 1 1; 0 .5; -1 1; -1 -1];
fh=@(p,x) 0.025 + 0.375*sum(p.^2,2);
[p,t]=distmesh2d(@dpoly,fh,0.025,[-1 -1; 1 1],pv,pv);
G = pebi(triangleGrid(p, t));
clf, plotGrid(G); axis equal tight off;
