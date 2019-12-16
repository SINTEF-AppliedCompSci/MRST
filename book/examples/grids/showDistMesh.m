%% Square domain with circular cutout, uniform mesh
% Use distmesh as included in the upr module. If you have your own distmesh
% module, use "mrstModule add distmesh" instead
mrstModule add upr;
linearize=false;
fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
[p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1],linearize);
G = triangleGrid(p, t);
clf, plotGrid(G); axis equal tight off;

%% Same domain, but graded grid
linearize=false;
fh=@(p) 0.05+0.3*dcircle(p,0,0,0.5);
[p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1],linearize);
G = triangleGrid(p, t);
clf, plotGrid(G); axis equal tight off;

%% Polyhedral domain, graded PEBI grid
linearize=false;
pv=[-1 -1; 0 -.5; 1 -1; 1 1; 0 .5; -1 1; -1 -1];
fh=@(p,x) 0.025 + 0.375*sum(p.^2,2);
[p,t]=distmesh2d(@dpoly,fh,0.025,[-1 -1; 1 1],pv,linearize,pv);
G = pebi(triangleGrid(p, t));
clf, plotGrid(G); axis equal tight off;
