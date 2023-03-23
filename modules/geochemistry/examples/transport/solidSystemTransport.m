close all;
clear;

%% load appropriate modules
mrstModule add ad-core ad-props ad-blackoil geochemistry
mrstVerbose off

%% Define the domain
% simple 1D cartesian domain 
G = cartGrid([10, 1, 1], [100, 1, 1]);
G = computeGeometry(G);
nc = G.cells.num;

%% Define the rock

rock.perm = 1*darcy*ones(G.cells.num, 1);
rock.poro = 0.3*ones(G.cells.num, 1);

%% Define the fluid
pRef = 0*barsa;
fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                           1000*kilogram/meter^3, 'c', 1e-10, 'cR', 4e-10, ...
                           'pRef', pRef);

%% Define the chemistry
% here we have the carbonate system with calcite precipitation

elements = {'O', 'H', 'Na*', 'Cl*', 'Ca*', 'C*'};

species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaOH', 'H2O*',...
             'Ca+2', 'CO3-2', 'HCO3-', 'CO2',...
             'CaCO3(s)'};

reactions ={'H2O  = H+  + OH- '       , 10^-14*mol/litre, ...
            'NaOH = Na+ + OH-'        , 10^10*mol/litre, ...
            'CaCO3(s) = CO3-2 + Ca+2' , 10^-8.48*(mol/litre)^2, ...
            'CO3-2 + H+ = HCO3-'      , 10^10.329/(mol/litre), ...
            'CO3-2 + 2*H+ = CO2 + H2O', 10^16.681/(mol/litre)};

chemsys = ChemicalSystem(elements, species, reactions);

% print the chemical model to the screen
chemsys.printChemicalSystem;

% instantiate chemical model
chemmodel = ChemicalModel(chemsys);

% initial chemistry
Nai = 1e-1;
Cli = Nai;
Cai = 1e-2;
Ci = 1e-3;
Hi = 1e-9;
H2Oi = 1;

inputConstraints = [Nai Cli Cai Ci Hi H2Oi]*mol/litre;
[initchemstate, initreport]= chemmodel.initState(repmat(inputConstraints, nc,1), 'charge', 'Cl');

% injected chemistry, we are injecting less Ca into the system
Naf = 1e-1;
Clf = Naf;
Caf = 1e-3;
Cf = 1e-3;
Hf = 1e-9;
H2Of = 1;

inputConstraints = [Naf Clf Caf Cf Hf H2Of]*mol/litre;
[injchemstate, injreport] = chemmodel.initState(inputConstraints, 'charge', 'Cl');

%% Define the initial state

initState = initchemstate;
initState.pressure = pRef*ones(nc,1);

%% instantiate the transport model
model = ChemicalTransportModel(G, rock, fluid, chemmodel);

%% Define the boundary conditions

% define the initial and injected fluid composition states
injfluidpart = injchemstate.elements;
initfluidpart = initchemstate.elements(end,:);

% sxtract the pore volumes from the domain
pv = poreVolume(G,rock);

% define source term at cell 1 for inflow boundary condition
src = [];
src = addSource(src, 1, pv(1)/day, 'sat', 1);

% give the fluid concentration at the inlet
src.elements    = injfluidpart;
src.logElements = log(injfluidpart);

% give dirchlet boundary condition at outlet
bc             = [];
bc             = pside(bc, G, 'east', 0*barsa, 'sat', 1);
bc.elements    = initfluidpart;        % (will not be used if outflow)
bc.logElements = log(initfluidpart);  % (will not be used if outflow)


%% Define the schedule
% it is recommended to ramp up the time steps

% ten time steps of 0.01 days followed by 100 steps of 1 day
schedule.step.val     = [0.01*day*ones(10, 1); 1*day*ones(10, 1);];
schedule.step.control = ones(numel(schedule.step.val), 1);
schedule.control      = struct('bc', bc, 'src', src, 'W', []);


%% Run the simulation

[~, states, scheduleReport] = simulateScheduleAD(initState, model, schedule);
mrstModule add mrst-gui

plotToolbar(G, states, 'field', 'elements:5', 'plot1d', true,'startplayback',true)
ylabel('Ca / mol/m^3')

%% Copyright notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2017 SINTEF ICT, Applied Mathematics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
