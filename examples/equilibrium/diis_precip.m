clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on


%% generate chemical system 

% define elements names
elements = {'Ba', 'Ca' 'O','S'};

% define chemical species
species = {'Ba+2+','Ca+2','SO4-2',...
            'BaSO4(s)','CaSO4(s)'};

% list chemical reactions         
reactions ={'BaSO4(s)  <-> Ba+2 + SO4-2 ',       1*mol/litre,...
            'CaSO4(s)  <-> Ca+2 + SO4-2',       0.67*mol/litre};
        


% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions);

chem.plotIter = false;

% print the chemical system
chem.printChemicalSystem;

%% solve the chemical system given inputs
n = 100;

C = 1e-3.*ones(n,1);
H2O = ones(n,1);
Alk = logspace(-5, -4,n)';

userInput = [C Alk]*mol/litre;

tic
[state, report, model] = chem.initState(userInput);
toc;

[state, chem] = chem.computeActivities(state);
[state, chem] = chem.computeChargeBalance(state);
[state, chem] = chem.computeSurfaceCharges(state);
[state, chem] = chem.computeSurfacePotentials(state);

plot(log10(Alk*litre/mol), log10(state.components*litre/mol));
xlabel('Alkalinity')
legend(chem.CompNames)