close all;
clear;

mrstModule add ad-core ad-props ad-blackoil geochemistry mrst-gui

mrstVerbose on

%% Define the grid

G = cartGrid([100, 1, 1], [10, 1, 1]);
G = computeGeometry(G);

plotGrid(G), view(3), axis tight

%% Define the rock

rock.perm = 1*darcy*ones(G.cells.num, 1);
rock.poro = 0.4*ones(G.cells.num, 1);

%% Define the fluid
pRef = 0*barsa;
tRef = 298*Kelvin;
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                           1000*kilogram/meter^3, 'c', 1e-10, 'cR', 4e-10, ...
                           'pRef', pRef);

%% Define the chemistry

%% generate chemical system 

% define elements names
elements = {'Ba*','Ca*','SO4*'};

% define species names
species = {'Ba+2','Ca+2','SO4-2',...
            'BaSO4(s)','CaSO4(s)'};
        

% list chemical reactions         
reactions ={'CaSO4(s)  <-> Ca+2 + SO4-2 ',       1*mol/litre,...
            'BaSO4(s)  <-> Ba+2 + SO4-2',        0.67*mol/litre};       

% list solid densities
solidDensities = {'CaSO4(s)', 3*mol/litre, 'BaSO4(s)',  2*mol/litre};

% instantiate the chemical model
chemModel = ChemicalModel(elements, species, reactions);

chemModel.plotIter = false;

% print the chemical system
chemModel.printChemicalSystem;

%% solve the chemical system given inputs
n =1;

Ba  = 0.5;
Ca  = 0.5;
SO4 = 1;
H = 10^-7;
H2O = 1;

injInput = [Ba Ca SO4]*mol/litre;

Ba  = 0.5;
Ca  = 0.5;
SO4 = 0.5;
H = 10^-4;
H2O = 1;

initInput = [Ba Ca SO4]*mol/litre;


% initial chemistry
[initchemstate, initreport]= chemModel.initState(initInput, 'solid', solidDensities);

% injected chemistry
[injchemstate, injreport] = chemModel.initState(injInput, 'solid', solidDensities);

%% Define the initial state

nc = G.cells.num;
fieldNames = fields(initchemstate);

for i = 1 : numel(fieldNames);
    initState.(fieldNames{i}) = repmat(initchemstate.(fieldNames{i}), nc, 1);
end

initState.pressure          = pRef*ones(nc,1);

% initState.solidDensities(1:2:end,:) = initState.solidDensities(1:2:end,:)*1.5;
%% Define the model

set(groot, 'defaultLineLineWidth', 3);
model = ChemicalTransportLogModel(G, rock, fluid, chemModel);
model.chemicalModel.nonlinearTolerance = 1e-12;
model.nonlinearTolerance = 1e-13;

fluidpart = model.fluidMat*((injchemstate.components)');
fluidpart = fluidpart';

model.plotIter = false;
%% Define the boundary conditions

src                  = [];
src                  = addSource(src, [1], [1/10].*meter^3/day, 'sat', 1);
src.masterComponents = fluidpart;
src.logmasterComponents = log(fluidpart);

bc                  = [];
bc                  = pside(bc, G, 'east', 0*barsa, 'sat', 1);
bc.masterComponents= [initchemstate.masterComponents];        % (will not used if outflow)
bc.logmasterComponents= [initchemstate.logMasterComponents];  % (will not used if outflow)


%% Define the schedule

schedule.step.val = [0.001*day*ones(4, 1); 0.1*day*ones(5, 1);1*day*ones(50, 1)];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('bc', bc, 'src', src, 'W', []);


%% Run the simulation

[~, states, scheduleReport] = simulateScheduleAD(initState, model, schedule);

% [ states ] = changeUnits( states, mol/litre );
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
