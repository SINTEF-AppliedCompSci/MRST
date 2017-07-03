close all;
clear;

mrstModule add ad-core ad-props ad-blackoil geochemistry mrst-gui
mrstVerbose on

%% Define the grid
nx = 100;
ny = 1;
nz = 20;

injInd = 2000-80;

G = cartGrid([nx, ny, nz], [100, 1, 20]);
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

elements = {'O', 'H', 'Na*','Cl*','Ca*','C'};

species = {'H+*', 'OH-','H2O*',...
           'Na+','Cl-','NaCl',...
           'Ca+2', 'CaCO3','CaHCO3+',...
           'CO3-2','HCO3-','CO2*',...
           '>SiO-', '>SiOH'};

reactions ={'H2O  <-> H+  + OH- ',      10^-14*mol/litre, ... 
            'NaCl <-> Na+ + Cl-',       10^1*mol/litre,...
            '>SiOH <-> >SiO- + H+',     10^-8*mol/litre,...
            'CO3-2 + 2*H+ <-> CO2 + H2O',   10^16.681/(mol/litre),...
            'HCO3- <-> H+ + CO3-2',      	10^-10.329*mol/litre,...
            'CaCO3 <-> Ca+2 + CO3-2',       10^3.224*mol/litre,...
            'Ca+2 + CO3-2 + H+ <-> CaHCO3+' 10^11.435/(mol/litre)^2};

geometry = [2*site/(nano*meter)^2 50e-3*meter^2/(gram) 5e3*gram/litre];
sioInfo = {geometry, 'tlm', [1 1e3], '>SiO-', [-1 0 0], '>SiOH', [0 0 0]};

surfaces = {'>SiO', sioInfo};

% instantiate chemical model
chemModel = ChemicalModel(elements, species, reactions, surfaces);
chemModel.printChemicalSystem;

%% initiate the chemical system
%'Na'    'Cl'    'Ca'    '>SiO'    'H+'    'H2O'    'CO2'

% initial chemistry
Na = 1e-3;
Cl = 1e-3;
Ca = 1e-4;
H = 1e-7;
H2O = 1;
CO2 = 1e-6;

inputConstraints = [Na, Cl, Ca, H, H2O, CO2]*mol/litre;
[initchemstate, initreport]= chemModel.initState(inputConstraints, 'charge', 'Cl');

% injected chemistry
Na = 1e-1;
Cl = 1e-1;
Ca = 1e-3;
H = 1e-4;
H2O = 1;
CO2 = 1e-3;

inputConstraints = [Na, Cl, Ca, H, H2O, CO2]*mol/litre;
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

initfluidpart = model.fluidMat*((initchemstate.components)');
initfluidpart = initfluidpart';

injfluidpart = model.fluidMat*((injchemstate.components)');
injfluidpart = injfluidpart';

model.plotIter = false;
%% Define the boundary conditions

src                  = [];
src                  = addSource(src, [(1:nx:nc)';injInd], [1/100*ones(nz,1);1/200].*meter^3/day, 'sat',ones(nz+1,1));
src.masterComponents = [repmat(initfluidpart,nz,1); injfluidpart];
src.logmasterComponents = log(src.masterComponents);

bc                  = [];
bc                  = pside(bc, G, 'east', 0*barsa, 'sat', 1);
bc.masterComponents= repmat(initfluidpart,nz,1);        % (will not used if outflow)
bc.logmasterComponents= log(bc.masterComponents);  % (will not used if outflow)


%% Define the schedule

schedule.step.val = [0.1*day*ones(5, 1);5*day*ones(10000, 1)];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control = struct('bc', bc, 'src', src, 'W', []);


%% Run the simulation

[~, states, scheduleReport] = simulateScheduleAD(initState, model, schedule);

save CO2Injector_run.mat -v7.3

[ states ] = changeUnits( states, mol/litre );

% 
v = VideoWriter('transport.avi');
open(v);

figure(1); box on;
xlabel('position')
ylabel('H+');
set(gca,'nextplot','replacechildren'); 
x = G.cells.centroids(:,1);
drawnow;
for i = 1 :10: numel(states)
    plot(x, -log10(states{i}.components(:,1)))
   frame = getframe(gcf);
   writeVideo(v,frame);
end

close(v);

x = 1:nx;
y = 1:nz;

pHmat = -log10(reshape(states{end}.components(:,1), nx,nz));

contourf(x,y, pHmat);
