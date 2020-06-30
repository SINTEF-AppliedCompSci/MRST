%% Mandel problem
%
% References : 
% @article{verruijt2013theory,
%          title={Theory and problems of poroelasticity},
%          author={Verruijt, Arnold},
%          journal={Delft University of Technology},
%          volume={71},
%          year={2013}
%         } (section 3.2)
%
% and
%
% @article{mikelic2014numerical,
%   title={Numerical convergence study of iterative coupling for coupled flow and geomechanics},
%   author={Mikeli{\'c}, Andro and Wang, Bin and Wheeler, Mary F},
%   journal={Computational Geosciences},
%   volume={18},
%   number={3-4},
%   pages={325--341},
%   year={2014},
%   publisher={Springer}
% }
%

clear all
% close all

%% Load required modules

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui mpsaw mpfa

%% Setup grid

% physdim = [20, 20] * meter;
physdim = [1, 1]*meter;
nx = 500; ny = 10;
resolution = [nx, ny];
G = cartGrid(resolution, physdim);
G = computeGeometry(G);

% Flow parameters
perm = 1;
muW  = 1;
poro = 1;

% fluid compressibility
cW = 0; 

% Bulk's modulus
K = 1;
% Poisson's ratio
nu = 0;
% Second Lamé parameter (also called Shear modulus and denoted G)
mu = 3/2*((1 - 2*nu)/(1 + nu))*K;
Gm = mu;

% First Lamé parameter
lambda = K - 2/3*Gm;
Gm = mu;

% Consolidation coefficient
cv = perm/muW*(K + 3*Gm); % see reference verruijt2013theory

% biot's coefficient
alpha = 1;

% force at top
top_force = 1;

% setup rock

rock.poro  = poro * ones(G.cells.num, 1);
rock.perm  = (perm/muW) * ones(G.cells.num, 1);
rock.alpha = alpha * ones(G.cells.num, 1);

% reference pressure on the side
pref = 0*barsa;

%% setup mechanics mech structure (with field prop and loadstruct)

lambda = lambda*ones(G.cells.num, 1);
mu = mu*ones(G.cells.num, 1);
mechprop = struct('lambda', lambda, 'mu', mu);

[tbls, mappings] = setupStandardTables(G);

% We recover the top, bottom and lateral faces using function pside
dummy = 0;
bc = pside([], G, 'Ymax', dummy); 
topfaces = bc.face;
bc = pside([], G, 'Xmin', dummy); 
leftfaces = bc.face;
bc = pside([], G, 'Xmax', dummy); 
rightfaces = bc.face;
bc = pside([], G, 'Ymin', dummy); 
bottomfaces = bc.face;

lateralfaces = leftfaces;

% At the bottom, we have rolling condition in x-direction
bottomlinform  = repmat([0, 1], numel(bottomfaces), 1);
% On the lateral walls, we have rolling condition in y-direction
laterallinform = repmat([1, 0], numel(lateralfaces), 1);
% At the top we have an unknown y-displacement but which is constant in x-direction. We incorporate it in the Dirichlet
% condition and add in the equations an extra variable, see mandelEquations
toplinform = repmat([0, 1], numel(topfaces), 1);

linform  = [bottomlinform; laterallinform; toplinform];
extfaces = [bottomfaces; lateralfaces; topfaces];
bcvals   = zeros(numel(extfaces), 1);

bc = struct('linform'    , linform , ...
            'extfaces'   , extfaces, ...
            'linformvals', bcvals);

bc = setupFaceBC(bc, G, tbls);

loadstruct.bc = bc;

nodefacecoltbl = tbls.nodefacecoltbl;
cellcoltbl = tbls.cellcoltbl;
loadstruct.extforce = zeros(nodefacecoltbl.num, 1); 
loadstruct.force = zeros(cellcoltbl.num, 1);

% Setup mech structure 
mech.prop = mechprop;
mech.loadstruct = loadstruct;

%% Setup flow parameters (with field c and bcstruct)

fluid.c = cW;
fluid.src = [];

% Setup boundary conditions for flow

bcfaces = rightfaces;
bcvals = pref*ones(numel(bcfaces));

nodefacetbl = tbls.nodefacetbl;

bcfacetbl.faces = bcfaces;
bcfacetbl = IndexArray(bcfacetbl);
bcnodefacetbl = crossIndexArray(nodefacetbl, bcfacetbl, {'faces'});

map = TensorMap();
map.fromTbl = bcfacetbl;
map.toTbl = bcnodefacetbl;
map.mergefds = {'faces'};
map = map.setup();

bcvals = map.eval(bcvals);

bcdirichlet = struct('bcnodefacetbl', bcnodefacetbl, ...
                     'bcvals', bcvals);
bcneumann = [];

bcstruct = struct('bcdirichlet', bcdirichlet, ...
                  'bcneumann'  , bcneumann);

fluid.bcstruct = bcstruct;

%% Setup Biot model

model = MandelModel(G, rock, fluid, mech, topfaces);
model = model.validateModel();

%% Setup schedule

tsteps = 100;
duration = 1;
t = duration/tsteps*ones(tsteps, 1);
tt = [1; 1 + cumsum(t)];
trepvals = [1e-5; 0.01; 0.1; 0.5; 1];
trep = 1 + trepvals/cv;
tt = [trep; tt];
tt = uniquetol(tt, 1e-9);
tt = [0; tt];
t = diff(tt);
schedule.step.val = t;

%%

schedule.step.control = 2*ones(numel(schedule.step.val), 1);
schedule.step.control(1) = 1;
schedule.control(1) = struct('W', [], 'avgtopforce', 0);
schedule.control(2) = struct('W', [], 'avgtopforce', top_force);

%% Setup initial state
clear initState;
% fluid
initState.pressure = zeros(G.cells.num, 1);
nlf = size(bcstruct.bcdirichlet.bcvals, 1);
initState.lambdafluid = zeros(nlf, 1);
% mech
cellcoltbl = tbls.cellcoltbl;
initState.u = zeros(cellcoltbl.num, 1);
nlm = size(loadstruct.bc.linformvals, 1);
initState.lambdamech = zeros(nlm, 1);
initState.avgtopforce = 0;
initState.vd = 0;

solver = NonLinearSolver('maxIterations', 100);
[wsol, states] = simulateScheduleAD(initState, model, schedule, 'nonlinearsolver', solver);

%% Mandel plot (pressure profile for some selected values)

figure
clf
ind = (1 : nx)';
xc = G.cells.centroids(ind, 1);
hold on

tt = cumsum(schedule.step.val);
legends = {};
for i = 1 : numel(states);
    [lia, locb] = ismembertol(tt(i), trep, 1e-8);
    if lia
        p = states{i}.pressure;
        p = p(ind);
        plot(xc, p);
        legends{end + 1} = sprintf('%g', trepvals(locb));
    end
end
legend(legends{:});

%% plotting

figure 
clf
plotToolbar(G, states);

figure
clf
ind = 1;
pmid = cellfun(@(state) state.pressure(ind), states);
tt = cumsum(schedule.step.val);
plot(tt, pmid, '*');
title('Pressure evolution at left edge')



