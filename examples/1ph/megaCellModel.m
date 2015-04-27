%% Simulation of a Mega-cell Model
% The purpose of this example is to show how one can setup MRST so that it
% is capable of solving models with multi-million cells. To this end, we
% will use the two-point flux-approximation scheme combined with a highly
% efficient algebraic multigrid method (AGMG) to solve the single-phase
% pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\frac{K}{\mu}\nabla p,$$
%
% for a flow driven by Dirichlet and Neumann boundary conditions. Our
% geological model will be simple a Cartesian grid with anisotropic,
% homogeneous permeability.


%% Define geometry
% Construct a Cartesian grid of size 200-by-100-by-100 (=2,000,000) cells
% in a box of dimensions 10-by-10-by-4 metres. We have been able to run
% this problem with a peak memory use of slightly less than 4GB.
nx = 200; ny = 100; nz = 100;
G = cartGrid([nx, ny, nz], [10, 10, 4]);


%% Process geometry
% The computation of geometric information (centroids and volumes of the
% cells and centroids, normals, and areas for the faces) is significantly
% faster if one uses the C-accelerated routine instead of the standard call
% to computeGeometry(G)
try
   require libgeometry incomp
catch
   mrstModule add libgeometry incomp
end
G = mcomputeGeometry(G);

%% Set rock and fluid data
% The only parameters in the single-phase pressure equation are the
% permeability $K$ and the fluid viscosity $\mu$.
rock.perm = repmat([1000, 100, 10].* milli*darcy(), [G.cells.num, 1]);
T         = computeTrans(G, rock);

gravity off;
fluid     = initSingleFluid('mu' ,    1*centi*poise     , ...
                            'rho', 1014*kilogram/meter^3);

%% Set boundary conditions
% Here, we impose Neumann conditions (flux of 1 m^3/day) on the global
% left-hand side. Similarly, we set Dirichlet boundary conditions p = 0 on
% the global right-hand side of the grid. For a single-phase
% flow, we need not specify the saturation at inflow boundaries and the
% fluid composition over outflow faces (here, right) is ignored by pside.
resSol = initResSol(G, 0.0);
bc     = fluxside([], G, 'LEFT',  1*meter^3/day());
bc     = pside   (bc, G, 'RIGHT', 0);


%% Solve the linear system
% Solve linear system construced from T and bc to obtain solution for flow
% and pressure in the reservoir. When working with large models, one cannot
% use the standard MLDIVIDE ('\') solver in MATLAB. Here, we use the AGMG
% algebraic multigrid solver.
mrstModule add agmg
resSol = incompTPFA(resSol, G, T, fluid, 'bc', bc, ...
                    'LinSolve', @(A,b) agmg(A,b,1));

clf
plotCellData(G, convertTo(resSol.pressure(1:G.cells.num), barsa()));
title('Cell Pressure [bar]')
xlabel('x'), ylabel('y'), zlabel('Depth');
view(3); camproj perspective; axis tight;
colorbar

%%
displayEndOfDemoMessage(mfilename)
