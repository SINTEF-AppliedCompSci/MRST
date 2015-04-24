if 1
    [grdecl, dataset, petroinfo] = getAtlasGrid('Utsirafm','coarsening',5);
    petrodata=petroinfo{1};

    G      = processGRDECL(grdecl{1});
else
    G = cartGrid([10,10,1], [10, 10, 1]*1000*meter);
end
G      = computeGeometry(G(1));
Gt     = topSurfaceGrid(G);

clear rock
rock.perm = 1*darcy*ones(G.cells.num, 1);
rock.poro = ones(G.cells.num, 1);
rho = [760 1200] .* kilogram/meter^3;

%%

[fluid, system, sol] = fluidFactoryVE(Gt, rock, 'fullyImplicit', true, 'dissolution', true);

solver = solverFactoryVE(Gt, fluid, rock, system, 'transportTolerance',1e-6);




res = trapAnalysis(Gt, true);
trees = maximizeTrapping(Gt, 'res', res);

[cell, largestVol, allFaces, point] = findOptimalInjectionPoint(Gt, res);
N = 5;
points = [];
for i = 1:N
    root = trees(i).root;
    cell = find(res.traps == root, 1);
    
    points = [points; ...
              Gt.cells.centroids(cell, :)];
end
[W, W_ve, bc] = wellSetup(Gt, rock, points, [100*mega*1000*kilogram]/N, rho, 'radius', .15);


mrstVerbose off
gravity z on 
% sol = initResSolVE(Gt, 300*barsa(), 0);
% sol.pressure= rho(2)*norm(gravity()).*Gt.cells.z;
% sol.wellSol = initWellSol(W, 300*barsa());
% sol.s = zeros(Gt.cells.num, 1);
% sol.rs=zeros(Gt.cells.num,1);
% sol.sGmax = sol.s;
% 
% for i = 1:10
%     i
%     [sol, its, conv] = solver(sol, W_ve, bc, 1*year);
% %     sol = solver(sol, [], bc, 1*year);
% end
%
close all
[states, timesteps] = solveMigrationGeneric(Gt, sol, W_ve, bc, fluid, solver,...
    'cutTimesteps', true, 'dt', 0.1*year);


%%
close all
plotToolbar(G, states)
axis tight off