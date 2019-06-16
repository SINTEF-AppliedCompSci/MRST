mrstModule add ad-core ad-props incomp mrst-gui mpfa postprocessing

%% Define and process geometry
% Construct a Cartesian grid of size 10-by-10-by-4 cells, where each cell
% has dimension 1-by-1-by-1. Because our flow solvers are applicable for
% general unstructured grids, the Cartesian grid is here represented using
% an unstructured formate in which cells, faces, nodes, etc. are given
% explicitly.
dimcase = 3;
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

clear vecs
caseno = 1;

% mpfa - jostein
T_mpfa = computeMultiPointTrans(G, rock, 'eta', eta);
state = initResSol(G, 0, 1);
state = incompMPFA(state, G, T_mpfa, fluid, 'bc', bc);
p = state.pressure;
vec = [z, p];
vecs{caseno} = sortrows(vec);
titles{caseno} = 'mpfa - jostein';
caseno = caseno + 1;

% mpfa - standard
mpfastruct = computeMultiPointTrans2(G, rock, 'eta', eta, 'verbose', true);
state = incompMPFA2(G, mpfastruct, 'bc', bc);
p = state.pressure;
vec = [z, p];
vecs{caseno} = sortrows(vec);
titles{caseno} = 'mpfa - standard';
caseno = caseno + 1;

% mpfa - block
mpfastruct = computeMultiPointTrans2(G, rock, 'eta', eta, 'blocksize', ...
                                      blocksize, 'verbose', true);
state = incompMPFA2(G, mpfastruct, 'bc', bc);
p = state.pressure;
vec = [z, p];
vecs{caseno} = sortrows(vec);
titles{caseno} = 'mpfa - block';
icaseno = caseno + 1;

% mpfa - new block
mpfastruct = blockComputeMultiPointTrans(G, rock, 'eta', eta, 'blocksize', ...
                                         blocksize, 'verbose', true);
state = incompMPFA2(G, mpfastruct, 'bc', bc);
p = state.pressure;
vec = [z, p];
vecs{caseno} = sortrows(vec);
titles{caseno} = 'mpfa - new block';
caseno = caseno + 1;


for i = 1 : numel(vecs)
    figure(i)
    clf
    plot(vecs{i}(:, 1), vecs{i}(:, 2));
    xlabel('z');
    ylabel('pressure');
    title(titles{i});
end


return

%% setup tracer bc

conc = ones(numel(bottomfaces), 1);
bctracer = setupBC(G, bottomfaces, conc);

%% Setup schedule
time  = 2.8e8*day;
% time  = 10*day;
dt    = 1e6*day;
dtvec = rampupTimesteps(time, dt, 20);
schedule = simpleSchedule(dtvec, 'bc', bctracer);

%% Setup tracer model

tracermodel = TracerModel(G, rock, fluid, 'tracerNames', {'tracer'});
state0.tracer = zeros(G.cells.num, 1);

[~, states] = simulateScheduleAD(state0, tracermodel, schedule);

%% plotting
clf, plotToolbar(G, states)