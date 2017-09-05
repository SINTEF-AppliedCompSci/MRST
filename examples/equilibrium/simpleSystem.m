%% create and solve chemical system
% clear;
close all;

% add autodiff core module
mrstModule add ad-core geochemistry
mrstVerbose off
%% generate chemical system object

elements = {'O', 'H', 'Na*', 'Cl*'};
species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaCl', 'H2O*'};
reactions = {'H2O  <-> H+  + OH- ', 10^-14*mol/litre, ...
             'NaCl <-> Na+ + Cl-',  10^1*mol/litre};

% instantiate chemical model
chem = ChemicalModel(elements, species, reactions);

chem.printChemicalSystem;

n = 500;
userInput = [1e-2.*ones(n,1) 1e-2.*ones(n,1) logspace(-5,-9, n)' 1.*ones(n,1)]*mol/litre;

[state, report, model] = chem.initState(userInput);

%% process

[state, chem] = chem.computeActivities(state);

% take out relevant values from state

state = changeUnits(state, {'masterComponents','components'}, mol/litre );

pH = -log10(getProp(chem, state, 'aH+'));

figure; hold on; box on;

for i = 1 : numel(chem.CompNames)
    v = getProps(chem, state, chem.CompNames{i});
    plot(pH, v)
    
end
xlabel('pH')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
legend(chem.CompNames)

