clear; 
close all;

%% Add geochemistry and automatic differentiation modules
mrstModule add ad-core geochemistry
mrstVerbose on

%% instantiate chemical model
% here we look at the pe dependent speciation of nitrogen. We are choosing
% total nitrogen, electron, sodium, H+ and H2O concentrations as inputs. 

elements = {'O', 'H', 'N*', 'e*','Na*'};

species = {'H+*', 'OH-', 'H2O*',...
            'NH4+', 'NO3-', 'e-','N2', 'NO2-','Na+'};

reactions ={'H2O  = H+  + OH-',                       10^-14*mol/litre,... 
            'NO3- + 10*H+ + 8*e- = NH4+ + 3*H2O',     10^119./(mol/litre)^15,...
            'NO3- + 2*H+ + 2*e- = NO2- + H2O',        10^28./(mol/litre)^3,...
            '2*NO3- + 12*H+ + 10*e- = N2 + 6*H2O',    10^(2*103)./(mol/litre)^17};
        
chemsys = ChemicalSystem(elements, species, reactions);

% print the chemical system
chemsys.printChemicalSystem;

% Setup model
chemmodel = ChemicalModel(chemsys);


%% define input values
% here we look at the speciation as a function of a range in electron
% concentration
n = 100;

N = 1e-3*ones(n,1);
e = logspace(-15, 10, n)';
H = 1e-7*ones(n,1);
H2O = ones(n,1);
Na = 1e-2*ones(n,1);

%% solve the chemical system
[state, report] = chemmodel.initState([N e Na H H2O]*mol/litre,'charge','H+');

%% calculate acitivities and charge balance
[state, chemmodel] = chemmodel.updateActivities(state);
[state, chemmodel] = chemmodel.updateChargeBalance(state);

state = changeUnits(state, {'species', 'activities'}, mol/litre);

pe = -log10(chemmodel.getProp(state, 'ae-'));

%% visualize the results
figure;
plot(pe, state.species, 'linewidth', 2);
set(gca, 'yscale', 'log');
ylabel('concentration [mol/litre]');
xlabel('pe');
legend(chemsys.speciesNames, 'location', 'southeast');

figure;
plot(pe, state.chargeBalance, 'linewidth', 2);
ylabel('charge balance ');
xlabel('pe');

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