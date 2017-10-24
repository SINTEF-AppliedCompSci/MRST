clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on


%% generate chemical system 

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
chem = ChemicalModel(elements, species, reactions, 'comb', combinations);

% print the chemical system
chem.printChemicalSystem;

%% solve the chemical system given inputs
n = 1000;

C = 1e-4.*ones(n,1);
H2O = ones(n,1);
Alk = linspace(-10^-3, 10^-3,n)';

userInput = [C H2O Alk]*mol/litre;

tic
[state, report, model] = chem.initState(userInput);
toc;

state = changeUnits(state, {'species'}, mol/litre);

figure; hold on;
plot(Alk, state.species, 'linewidth', 2);
set(gca,'yscale', 'log')
xlabel('Alkalinity [mol/litre]');
ylabel('concentration [mol/litre]');

legend(chem.speciesNames)