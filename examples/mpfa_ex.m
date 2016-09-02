% Example setting up and solving a simple 2D flow problem using MPFA. 
%
% For information on MRST-functions, confer the MRST documentation at 
%   http://www.sintef.no/projectweb/mrst/

clear
% Initialize a Cartesian grid.
% NOTE: If this does not work, it is most probably because MRST is not in
% your path.
G = computeGeometry(cartGrid([2, 3], [1, 1]));
% Plot the grid using MRST functions
figure; plotGrid(G)

% Initialize permeability. Syntax and storage as a structure is chosen for
% compatibility with the equivalent method in MRST.
rock.perm = ones(G.cells.num, 1) * [1, 3, .5];% Perm: (kx=1, ky=3, kxy=.5)

% Indicies of boundary faces are hard-coded according to the current grid
% size. MRST also has some functionality for picking out boundary faces, or
% we could have searched based on coordinates
xmin_faces = [1, 4, 7]; % Indicies of faces on x=0
xmax_faces = [3, 6, 9]; % Indicies of faces at x=1

% Boundary conditions: Dirichlet on x=0, the rest will default to Neumann
% Note that the value given for the BC is somehow irrelevant (but needed 
% for compatibility with MRST). We will instead assign the boundary values
% directly to bc_vals, below.
bc = addBC([], xmin_faces, 'pressure', 0); 

% The parameter invertBlocks describes which function is used for inverting
% local systems. This can be either 'matlab' or 'mex'. The mex option can
% give substantial speedup for large systems, in particular for elasticity,
% but it can also be memory-demanding.
mpfa_discr = mpfa(G, rock, [], 'invertBlocks', 'matlab', 'bc', bc);

% Assign boundary conditions
bc_vals = zeros(G.faces.num, 1);
bc_vals(xmin_faces) = 1:3;
bc_vals(xmax_faces) = -0.2;

% Flow is driven by boundary conditions
rhs = -mpfa_discr.div * mpfa_discr.boundFlux * bc_vals;
pressure = mpfa_discr.A \ rhs;
% Plot pressure
plotCellData(G, pressure)
% Flux consists of contributions from pressure and bcs
flux = mpfa_discr.F * pressure + mpfa_discr.boundFlux * bc_vals;