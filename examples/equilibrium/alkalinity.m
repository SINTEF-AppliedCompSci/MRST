clear;
close all;

%% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on

%% define the chemical model
% here we choose CO3, H2O and alkalaninity as inputs

% define elements names
elements = {'H', 'CO3*','O'};

% define chemical species
species = {'H+','H2CO3','HCO3-','CO3-2','H2O*','OH-'};

% list chemical reactions         
reactions ={'H2CO3  = H+ + HCO3- ',         10^2*mol/litre,...
            'HCO3- = H+ + CO3-2',           10^-1*mol/litre,...
            'H2O = H+ + OH-',               10^-14*mol/litre};
        
% define any linear combinations
combinations = {'Alk*', 'HCO3- + 2*CO3-2 + OH- - H+' };

% instantiate the chemical model
chemsys = ChemicalSystem(elements, species, reactions, 'comb', combinations);

% print the chemical system
chemsys.printChemicalSystem;

% Setup model
chemmodel = ChemicalModel(chemsys);


%% solve the chemical system given inputs
% here we vary alkalinity and keep carbon and water content constant
n = 1000;

C = 1e-4.*ones(n,1);
H2O = ones(n,1);
Alk = linspace(-10^-3, 10^-3,n)';

userInput = [C H2O Alk]*mol/litre;

tic
[state, report, model] = chemmodel.initState(userInput);
toc;

state = changeUnits(state, {'species'}, mol/litre);

% plot the results
figure; hold on;
plot(Alk, state.species, 'linewidth', 2);
set(gca,'yscale', 'log')
xlabel('Alkalinity [mol/litre]');
ylabel('concentration [mol/litre]');

legend(chemsys.speciesNames)

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
