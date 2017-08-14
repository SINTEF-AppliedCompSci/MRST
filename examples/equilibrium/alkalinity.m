clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on


%% generate chemical system 

% define elements names
elements = {'H', 'CO3*'};

% define chemical species
species = {'H+','H2CO3','HCO3-','CO3-2'};

% list chemical reactions         
reactions ={'H2CO3  <-> H+ + HCO3- ',       10^2*mol/litre,...
            'HCO3- <-> H+ + CO3-2',         10^-1*mol/litre};
        
% define any linear combinations
combinations = {'Alk*', 'HCO3- + 2*H2CO3' };

% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions, 'comb', combinations);

chem.plotIter = false;

% print the chemical system
chem.printChemicalSystem;

%% solve the chemical system given inputs
n = 100;

C = 1e-3.*ones(n,1);
H2O = ones(n,1);
Alk = logspace(-10, -4,n)';

userInput = [C Alk]*mol/litre;

tic
[state, report, model] = chem.initState(userInput);
toc;

[state, chem] = chem.computeActivities(state);
[state, chem] = chem.computeChargeBalance(state);
[state, chem] = chem.computeSurfaceCharges(state);
[state, chem] = chem.computeSurfacePotentials(state);

plot(log10(Alk), log10(state.components*litre/mol));
xlabel('Alkalinity')
legend(chem.CompNames)