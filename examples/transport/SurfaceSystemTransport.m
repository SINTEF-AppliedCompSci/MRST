close all;
clear;

mrstModule add ad-core ad-props ad-blackoil geochemistry mrst-gui


%% Define the grid

G = cartGrid([100, 1, 1], [10, 1, 1]);
G = computeGeometry(G);

plotGrid(G), view(3), axis tight

%% Define the rock

rock.perm = 1*darcy*ones(G.cells.num, 1);
rock.poro = ones(G.cells.num, 1);

%% Define the fluid
pRef = 0*barsa;
tRef = 298*Kelvin;
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                           1000*kilogram/meter^3, 'c', 1e-10, 'cR', 4e-10, ...
                           'pRef', pRef);

%% Define the chemistry

elements = {'O', 'H', 'Na*','Cl*'};

species = {'H+*', 'OH-', 'Na+', 'H2O*', '>SiO-', '>SiOH','NaCl','Cl-'};

reactions ={'H2O  <-> H+  + OH- ',      10^-14*mol/litre, ... 
            'NaCl <-> Na+ + Cl-',       10^1*mol/litre,...
            '>SiOH <-> >SiO- + H+',     10^-8*mol/litre};

geometry = [2*site/(nano*meter)^2 50e-3*meter^2/(gram) 5e3*gram/litre];
sioInfo = {geometry, 'tlm', [1 1e3], '>SiO-', [-1 0 0], '>SiOH', [0 0 0]};

surfaces = {'>SiO', sioInfo};

% instantiate chemical model
chemModel = ChemicalModel(elements, species, reactions, surfaces);
chemModel.printChemicalSystem;

% initial chemistry

inputConstraints = [1e-3 1e-3 1e-5 1]*mol/litre;
[initchemstate, initreport]= chemModel.initState(inputConstraints, 'charge', 'Cl');

% injected chemistry
inputConstraints = [1e-1 1e-1 1e-9 1]*mol/litre;
[injchemstate, injreport] = chemModel.initState(inputConstraints, 'charge', 'Cl');

%% Define the initial state

nc = G.cells.num;
initState.components          = repmat(initchemstate.components, nc, 1);
initState.masterComponents    = repmat(initchemstate.masterComponents, nc, 1);
initState.logcomponents       = repmat(initchemstate.logcomponents, nc, 1);
initState.logmasterComponents = repmat(initchemstate.logmasterComponents, nc, 1);

initState.pressure          = pRef*ones(nc,1);

%% Define the model

set(groot, 'defaultLineLineWidth', 3);
model = ChemicalTransportLogModel(G, rock, fluid, chemModel);
model.chemicalModel.nonlinearTolerance = 1e-12;
model.nonlinearTolerance = 1e-13;

fluidpart = model.fluidMat*((injchemstate.components)');
fluidpart = fluidpart';

model.plotIter = true;
%% Define the boundary conditions

src                  = [];
src                  = addSource(src, [1], [1/10].*meter^3/day, 'sat', 1);
src.masterComponents = fluidpart;
src.logmasterComponents = log(fluidpart);

bc                  = [];
bc                  = pside(bc, G, 'east', 0*barsa, 'sat', 1);
bc.masterComponents= [initchemstate.masterComponents];        % (will not used if outflow)
bc.logmasterComponents= [initchemstate.logmasterComponents];  % (will not used if outflow)


%% Define the schedule

schedule.step.val = [0.01*day*ones(200, 1); 1*day*ones(100, 1)];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('bc', bc, 'src', src, 'W', []);


%% Run the simulation

[~, states, scheduleReport] = simulateScheduleAD(initState, model, schedule);

[ states ] = changeUnits( states, mol/litre );
% 
% plotToolbar(G, states,'plot1d', true, 'log10', true);
% 
% v = VideoWriter('transport.avi');
% open(v);
% 
% figure(1); box on;
% xlabel('position')
% ylabel('pH');
% set(gca,'nextplot','replacechildren'); 
% x = G.cells.centroids(:,1);
% ylim([8.9 10]);
% drawnow;
% for i = 1 : numel(states)
%     plot(x, -log10(states{i}.components(:,1)))
%    frame = getframe(gcf);
%    writeVideo(v,frame);
% end
% 
% close(v);
