clear;
close all;

%% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on

%% instantiate the chemical system 
% here we look at carbonate speciation with calcite and CO2(g) as non
% aquoeus phases. We choose the partial pressure of CO2 as an input by
% appending an asterisk to the name

elements = {'O', 'H', 'Na*', 'Cl*', 'Ca*', 'C'};

species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaOH', 'H2O*',...
             'Ca+2', 'CO3-2', 'HCO3-', 'CO2',...
             'CaCO3(s)', 'CO2(g)*'};

reactions ={'H2O  = H+  + OH- ',            10^-14*mol/litre, ...
            'NaOH = Na+ + OH-',             10^10*mol/litre,...
            'CaCO3(s) = CO3-2 + Ca+2',      10^-8.48*(mol/litre)^2,...
            'CO3-2 + H+ = HCO3-',           10^10.329/(mol/litre),...
            'CO3-2 + 2*H+ = CO2 + H2O',     10^16.681/(mol/litre),...
            'CO2 = CO2(g)',                 10^1.468*atm/(mol/litre)};

chemsys = ChemicalSystem(elements, species, reactions);

% print the chemical system
chemsys.printChemicalSystem;

% Setup model
chemmodel = ChemicalModel(chemsys);


%% specify inputs
% here we vary the partial pressure of CO2 from 1e-3 atm to 3 atm
n = 100;

Na = 1e-1*ones(n,1)*mol/litre;
Cl = 1e-1*ones(n,1)*mol/litre;
H = 1e-4*ones(n,1)*mol/litre;
H2O = ones(n,1)*mol/litre;
Ca = 1e-2*ones(n,1)*mol/litre;
CO2 = logspace(-3, 1, n)'*atm;

%% solve chemical system given inputs
state = chemmodel.initState([Na, Cl, Ca, H, H2O, CO2], 'charge','H+');

%% process data
state = changeUnits(state, {'elements','species','partialPressures'}, [mol/litre, mol/litre, atm] );

[state, chemmodel] = chemmodel.updateActivities(state);

%% plot the results
figure; hold on; box on;
plot(CO2/atm, state.species, 'linewidth', 2)
set(gca, 'xscale','log');

xlabel('CO_2(g) [atm]')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
legend(chemsys.speciesNames)

figure; hold on; box on;
plot(CO2/atm, state.saturationIndicies, 'linewidth', 2)
set(gca, 'xscale','log');

xlabel('CO_2(g) [atm]')
ylabel('saturation indicies [-]');
set(gca, 'yscale', 'log');
legend(chemsys.solidNames)

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
