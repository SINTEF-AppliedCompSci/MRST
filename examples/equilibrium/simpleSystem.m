%% create and solve chemical system
% clear;
close all;

mrstModule add geochemistry
mrstVerbose on
%% generate chemical system object

elements = {'O', 'H', 'Na*', 'Cl*'};

species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaCl', 'H2O*'};

reactions = {'H2O  = H+  + OH- ', 10^-14*mol/litre, ...
             'NaCl = Na+ + Cl-',  10^1*mol/litre};

chem = ChemicalModel(elements, species, reactions); %, 'combinations', comboComps);

chem.printChemicalSystem;
chem.plotIterations = true;
chem.nonlinearTolerance = 1e-14;

n = 100;

Na  = ones(n,1)*milli*mol/litre;
Cl  = ones(n,1)*milli*mol/litre;
H2O = ones(n,1)*mol/litre;
H   = logspace(-6, -8,n)'*mol/litre;

inputs = [Na, Cl, H, H2O];

[state, report, chem] = chem.initState(inputs, 'chargeBalance', 'Na');

%% process
state = changeUnits(state, {'elements','species'}, mol/litre );

[state, chem] = chem.computeActivities(state);

pH = -log10(getProp(chem, state, 'aH+'));

figure; hold on; box on;
plot(pH, state.species)

xlabel('pH')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
legend(chem.speciesNames)

