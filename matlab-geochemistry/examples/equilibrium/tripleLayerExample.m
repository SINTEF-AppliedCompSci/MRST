clear;
close all;

%% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on

%% generate chemical system 
% here we look at an amphoteric surface with the sorption of Na and Cl

% define elements names
elements = {'O', 'H', 'Na*', 'Cl*'};

% define chemical species
species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaCl', 'H2O*',...
             '>SO-', '>SOH', '>SOH2+', '>SONa', '>SOH2Cl'};

% list chemical reactions         
reactions ={'H2O  = H+  + OH- ',          10^-14*mol/litre, ...
            'NaCl = Na+ + Cl-',           10^1*mol/litre,...
            '>SOH = >SO- + H+',           10^-7.5*mol/litre,...
            '>SOH + H+ = >SOH2+',         10^2/(mol/litre),...
            '>SO- + Na+ = >SONa',         10^-1.9/(mol/litre),...
            '>SOH2+ + Cl- = >SOH2Cl',     10^-1.9/(mol/litre)};
        
% define the surface, here we choose a triple layer model to describe the
% surface. Outersphere complexes like >SONa and >SOH2Cl are defined by
% defining their charge contributions to the inner and outer helmholtz plane

geometry = [2*site/(nano*meter)^2 50e-3*meter^2/gram 5e3*gram/litre];
sioInfo = {geometry, 'tlm', [1 0.2]/meter^2,      '>SONa',   [-1 1], '>SOH2Cl',[1 -1]};
surfaces ={ '>SO', sioInfo };

% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions, 'surf', surfaces);

% print the chemical system
chem.printChemicalSystem;

%% define the inputs
% here we vary pH from 4 to 10
n = 500;

Na = 1e-2.*ones(n,1);
Cl = 1e-2*ones(n,1);
H2O = ones(n,1);
H = logspace(-4, -10, n)';

userInput = [Na Cl H H2O]*mol/litre;


[state, report, model] = chem.initState(userInput, 'chargeBalance', 'Na');

%% post processing
% to add the surface potentials and charges to the state structure we must
% call on functions within ChemicalModel
[state, chem] = chem.computeActivities(state);
[state, chem] = chem.computeChargeBalance(state);
[state, chem] = chem.computeSurfaceCharges(state);
[state, chem] = chem.computeSurfacePotentials(state);
[state, chem] = chem.computeAqueousConcentrations(state);
[state, chem] = chem.computeSurfaceConcentrations(state);

state = changeUnits(state, {'elements', 'species', 'activities'}, mol/litre);

pH = -log10(getProp(chem, state, 'aH+'));

%% plot the results

% plot the charge balance of the system
figure; hold on; box on;
plot(pH, state.chargeBalance,'-k')
xlabel('pH')
ylabel('charge [% of total ion concentration]');

% plot the concentration of elements on the surface
figure; hold on; box on;
plot(pH, state.surfaceConcentrations, 'linewidth',2)
xlabel('pH')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
legend(chem.surfaceConcentrationNames);

% plot the surface potentials, there are three surfaces in the triple layer
% model
figure; hold on; box on;
plot(pH, state.surfacePotentials, 'linewidth',2)
xlabel('pH')
ylabel('potential [volts]');
legend(chem.surfacePotentialNames);

% plot the surface charges, they will sum to zero
figure; hold on; box on;
plot(pH, state.surfaceCharges, 'linewidth',2)
xlabel('pH')
ylabel('charge density [C/m^2]');
legend(chem.surfaceChargeNames);

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