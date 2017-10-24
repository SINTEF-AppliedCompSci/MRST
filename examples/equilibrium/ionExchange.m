%% create and solve chemical system
clear;
close all;

mrstModule add geochemistry ad-core
mrstVerbose on
%% generate chemical system object

elements = {'O', 'H', 'Na*', 'Cl*'};

species = {'H+*', 'OH-', 'H2O*', '>YH', '>YNa', 'Na+', 'Cl-','NaOH'};

reactions = {'H2O  = H+  + OH- ',           10^-14*mol/litre, ...
             '>YH + Na+ = >YNa + H+',       10^10,...
             'NaOH = Na+ + OH-',            10^10*mol/litre};

surfInfo = {'>Y', {[1*site/(nano*meter)^2 1*meter^2/gram 1*gram/litre], 'ie'}};

chem = ChemicalModel(elements, species, reactions, 'surfaces', surfInfo);

%%
n = 100;

H2O = ones(n,1);
H   = logspace(-3, -11,n)';
Na  = 1e-2*ones(n,1);
Cl  = Na;

inputs = [Na, Cl, H, H2O]*mol/litre;

[state, report, chem] = chem.initState(inputs, 'chargeBalance', 'Na');

%% process
state = changeUnits(state, {'elements','species'}, mol/litre );

[state, chem] = chem.computeActivities(state);

pH = -log10(getProp(chem, state, 'aH+'));

figure; hold on; box on;
plot(pH, state.species, 'linewidth', 2)

xlabel('pH')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
legend(chem.speciesNames)
xlim([3 11])