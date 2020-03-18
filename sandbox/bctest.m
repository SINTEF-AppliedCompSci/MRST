clear all
close all

mrstModule add ad-core ad-props incomp mrst-gui postprocessing ...
    nfvm

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
rock = makeRock(G, 1e-3*darcy, 1);

fluid = initSingleFluid('mu' , 1, ...
                        'rho', 1);
Wtp = [];

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
bc2 = struct('face', [], 'type', {{}}, 'value', [], 'sat', []);
bc2.face = [bottomfaces; topfaces];
bc2.type = repmat({'pressure'}, [1, numel(bc2.face)]);
value = zeros(nbf + ntf, 1);
value(1 : nbf) = 1*barsa;
bc2.value = value;
bc = convertBC2FlowNTPFA(G, bc2);

%% Pressure run
titles = {};
z = G.cells.centroids(:, dim);
eta = 1/3;
blocksize = 10;

clear vecs fluxes
caseno = 1;

%% TPFA (also used for initial guess to NTPFA)
state0 = initResSol(G, 0, 1);
T = computeTrans(G, rock);
state = incompTPFA(state0, G, T, fluid, 'bc', bc2);
p0 = state.pressure;
%p0 = barsa*ones(G.cells.num, 1);
%rng(1)
%p0 = barsa*rand(G.cells.num, 1);
p              = state.pressure;
vec            = [z, p];
vecs{caseno}   = sortrows(vec);
fluxes{caseno} = state.flux;
titles{caseno} = 'tpfa';
caseno         = caseno + 1;



%% ntpfa kobaisi
tol = 1e-15;
mrstVerbose on
interpFace=findHAP(G,rock,bc);
OSflux=findOSflux(G,rock,bc,interpFace);
state=FlowNTPFA(G,bc,fluid,Wtp,OSflux,p0,tol,1000);
mrstVerbose off
p              = state.pressure;
vec            = [z, p];
vecs{caseno}   = sortrows(vec);
fluxes{caseno} = state.flux;
titles{caseno} = 'ntpfa-kobaisi';
caseno         = caseno + 1;

close all
for i = 1 : numel(vecs)
    figure
    plot(vecs{i}(:, 1), vecs{i}(:, 2), '.-');
    xlabel('z');
    ylabel('pressure');
    title(titles{i});
end

%% check flux computations by computing mass directly mass conservation in
%% each cell

nc = G.cells.num;
nf = G.faces.num;

N = G.faces.neighbors;

cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos));
cellfacetbl.faces = G.cells.faces(:, 1);
cellfacetbl.num   = numel(cellfacetbl.cells);

facetbl.faces = (1 : nf)';
facetbl.num   = nf;

sgn = zeros(cellfacetbl.num, 1);

faces = find(N(facetbl.faces, 1) > 0);
clear poscellfacetbl
poscellfacetbl.cells = N(faces, 1);
poscellfacetbl.faces = faces;
poscellfacetbl.num   = numel(poscellfacetbl.cells);
mappos = setupTableMapping(poscellfacetbl, cellfacetbl, {'cells', 'faces'});
sgn = sgn + mappos*ones(poscellfacetbl.num, 1);

faces = find(N(facetbl.faces, 2) > 0);
clear negcellfacetbl
negcellfacetbl.cells = N(faces, 2);
negcellfacetbl.faces = faces;
negcellfacetbl.num   = numel(negcellfacetbl.cells);
mapneg = setupTableMapping(negcellfacetbl, cellfacetbl, {'cells', 'faces'});
sgn = sgn - mapneg*ones(negcellfacetbl.num, 1);

tbl = cellfacetbl; % alias
div = sparse(tbl.cells, tbl.faces, sgn, nc, nf);

for i = 1 : numel(fluxes)
    figure
    q = div*(fluxes{i});
    plotCellData(G, q);
    title(['mass conservation - ', titles{i}]);
    colorbar
end
