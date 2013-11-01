[grdecl, dataset, petroinfo] = getAtlasGrid('Utsirafm','coarsening',5);
petrodata=petroinfo{1};
G      = processGRDECL(grdecl{1});
G      = computeGeometry(G(1));
Gt     = topSurfaceGrid(G);

clear rock
rock.perm = ones(G.cells.num, 1);
rock.poro = ones(G.cells.num, 1);
rho = [760 1200] .* kilogram/meter^3;

%%

[fluid, system] = fluidFactoryVE(Gt, rock);
solver = solverFactoryVE(Gt, fluid, rock, system);
[W, W_ve, bc] = wellSetup(Gt, rock, [1; 2], 1e9, rho);
 
%%
sol = initResSolVE(Gt, 300*barsa(), 0);
sol.pressure= rho(2)*norm(gravity()).*Gt.cells.z;
sol.wellSol = initWellSol(W, 300*barsa());
sol.s = zeros(Gt.cells.num, 1);
sol.rs=zeros(Gt.cells.num,1);
sol.sGmax = sol.s;

for i = 1:10
    sol = solver(sol, W_ve, bc, 1*year);
end