close all;
clear;

mrstModule add ad-core ad-props ad-blackoil geochemistry mrst-gui

mrstVerbose on

%% Define the grid

G = cartGrid([100, 1, 1], [10, 1, 1]);
G = computeGeometry(G);
nc = G.cells.num;
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

species = {'H+*', 'OH-', 'Na+', 'H2O*', 'NaCl','Cl-'};

reactions ={'H2O  = H+  + OH- ',      10^-14*mol/litre, ... 
            'NaCl = Na+ + Cl-',       10^1*mol/litre};


% instantiate chemical model
chemModel = ChemicalModel(elements, species, reactions);
chemModel.printChemicalSystem;

% initial chemistry

inputConstraints = [1e-1 1e-1 1e-5 1]*mol/litre;
[initchemstate, initreport]= chemModel.initState(repmat(inputConstraints, nc/2,1), 'charge', 'Cl');

% injected chemistry
inputConstraints = [1e-3 1e-3 1e-9 1]*mol/litre;
[injchemstate, injreport]= chemModel.initState(repmat(inputConstraints, nc/2,1), 'charge', 'Cl');

%% Define the initial state


initState = initchemstate;
initState.pressure          = pRef*ones(nc,1);

%% Define the model

set(groot, 'defaultLineLineWidth', 3);
model = ChemicalTransportLogModel(G, rock, fluid, chemModel);

injFluidPart = (model.fluidMat*injchemstate.species(end,:)')';

initFluidPart = (model.fluidMat*initchemstate.species(end,:)')';

model.plotIter = false;
%% Define the boundary conditions

src                  = [];
src                  = addSource(src, [1], [1/10].*meter^3/day, 'sat', 1);
src.masterComponents = injFluidPart;
src.logmasterComponents = log(injFluidPart);

bc                  = [];
bc                  = pside(bc, G, 'east', 0*barsa, 'sat', 1);
bc.masterComponents= initFluidPart;        % (will not used if outflow)
bc.logMasterComponents= log(initFluidPart);  % (will not used if outflow)


%% Define the schedule

schedule.step.val = [0.1*day*ones(5, 1); 1*day*ones(5, 1); 3*day*ones(50, 1)];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('bc', bc, 'src', src, 'W', []);


%% Run the simulation

[~, states, scheduleReport] = simulateScheduleAD(initState, model, schedule);

[ states ] = changeUnits( states, {'components', 'masterComponents'}, mol/litre );
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
