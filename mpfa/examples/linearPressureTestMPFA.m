%% Linear pressure test
%
%
% The MPFA method is exact for linear pressure field. We impose boundary
% condition such that the exact solution is linear. The grid is a twisted grid
% made from a Cartesian grid.
%
% We check the three implementations (Legacy, standard, block assembly)


mrstModule add ad-core ad-props incomp mrst-gui mpfa postprocessing

clear all

%% Define and process geometry
% Construct a Cartesian grid of size 10-by-10-by-4 cells, where each cell
% has dimension 1-by-1-by-1. Because our flow solvers are applicable for
% general unstructured grids, the Cartesian grid is here represented using
% an unstructured formate in which cells, faces, nodes, etc. are given
% explicitly.
dimcase = 2;
switch dimcase
  case 2
    nx = 10; ny = 10;
    G = cartGrid([nx, ny]);
  case 3
    nx = 4; ny = 3; nz = 5;
    G = cartGrid([nx, ny, nz]);
end
G = twister(G, 0.1); 
G = computeGeometry(G);
nc = G.cells.num;

%% Set rock and fluid data
% The only parameters in the single-phase pressure equation are the
% permeability $K$, which here is homogeneous, isotropic and equal 100 mD.
% The fluid has density 1000 kg/m^3 and viscosity 1 cP.
% We make a non diagonal rock tensor
% rock = makeRock(G, 1e-3*darcy, 1);
rock = makeRock(G, 1, 1);

fluid = initSingleFluid('mu' , 1, ...
                        'rho', 1);

gravity off

%% setup bc for pressure
dim = G.griddim;
zface = G.faces.centroids(:, dim);
maxz = max(zface);
minz = min(zface);
bottomfaces = find(G.faces.centroids(:, dim) > (maxz - 1e-4));
topfaces = find(G.faces.centroids(:, dim) < (minz + 1e-4));

nbf = numel(bottomfaces);
ntf = numel(topfaces);
bc = struct('face', [], 'type', {{}}, 'value', [], 'sat', []);
bc.face = [bottomfaces; topfaces];
bc.type = repmat({'pressure'}, [1, numel(bc.face)]);
value = zeros(nbf + ntf, 1);
value(1 : nbf) = 1*barsa;
bc.value = value;


%% Pressure run
titles = {};
z = G.cells.centroids(:, dim);
eta = 1/3;
blocksize = 10;

clear vecs fluxes
caseno = 1;

% mpfa - Legacy implementation
T_mpfa = computeMultiPointTransLegacy(G, rock, 'eta', eta);
state = initResSol(G, 0, 1);
state = incompMPFAlegacy(state, G, T_mpfa, fluid, 'bc', bc);
p              = state.pressure;
vec            = [z, p];
vecs{caseno}   = sortrows(vec);
fluxes{caseno} = state.flux;
titles{caseno} = 'mpfa - legacy';
caseno         = caseno + 1;

% mpfa - Standard
mpfastruct = computeMultiPointTrans(G, rock, 'eta', eta, 'verbose', true);
state = incompMPFAbc(G, mpfastruct, bc, 'outputFlux', true);
mpfastruct1 = mpfastruct;
p              = state.pressure;
vec            = [z, p];
vecs{caseno}   = sortrows(vec);
fluxes{caseno} = state.flux;
titles{caseno} = 'mpfa - standard';
caseno         = caseno + 1;

% mpfa - block assembly (necessary for large systems)
mpfastruct = blockComputeMultiPointTrans(G, rock, 'eta', eta, 'blocksize', ...
                                         blocksize, 'verbose', true);
mpfastruct2 = mpfastruct;
state = incompMPFAbc(G, mpfastruct, bc, 'outputFlux', true);
p              = state.pressure;
vec            = [z, p];
vecs{caseno}   = sortrows(vec);
fluxes{caseno} = state.flux;
titles{caseno} = 'mpfa - block';
caseno         = caseno + 1;

close all
for i = 1 : numel(vecs)
    figure
    plot(vecs{i}(:, 1), vecs{i}(:, 2), '*-');
    xlabel('z');
    ylabel('pressure');
    title(titles{i});
end
