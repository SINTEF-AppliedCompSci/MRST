mrstModule add ad-core ad-props incomp mrst-gui mpfa postprocessing eni

%% Define and process geometry
% Construct a Cartesian grid of size 10-by-10-by-4 cells, where each cell
% has dimension 1-by-1-by-1. Because our flow solvers are applicable for
% general unstructured grids, the Cartesian grid is here represented using
% an unstructured formate in which cells, faces, nodes, etc. are given
% explicitly.
nx = 4; ny = 3; nz = 5;
G = cartGrid([nx, ny, nz]);
% G = twister(G, 0.1); 
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

clear states
mpfastruct = computeMultiPointTrans2(G, rock, 'eta', 1/3);
states{1} = incompMPFA2(G, mpfastruct, 'bc', bc);

p = states{1}.pressure;
z = G.cells.centroids(:, dim);
vec = [z, p];
vec = sortrows(vec);

figure
plot(vec(:, 1), vec(:, 2));
xlabel('z');
ylabel('pressure');

T_mpfa = computeMultiPointTrans(G, rock,'eta',1/3);
states{2} = initResSol(G, 0, 1);
states{2} = incompMPFA(states{2}, G, T_mpfa, fluid, 'bc', bc);

p = states{2}.pressure;
z = G.cells.centroids(:, dim);
vec = [z, p];
vec = sortrows(vec);

figure
plot(vec(:, 1), vec(:, 2));
xlabel('z');
ylabel('pressure');

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