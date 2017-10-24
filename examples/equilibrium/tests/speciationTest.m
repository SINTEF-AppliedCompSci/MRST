%% create and solve chemical system
clear;
close all;

mrstModule add geochemistry ad-core
mrstVerbose on
%% generate chemical system object

elements = {'O', 'H', 'Na*', 'Cl*', 'Ca*'};

species = {'H+*', 'OH-', 'H2O*',...
           'Na+', 'Cl-', 'NaOH',...
           'Ca+2', 'CaOH+'};

reactions = {'H2O  = H+  + OH- ',       10^-14*mol/litre, ...
             'NaOH = Na+ + OH-',        10^10*mol/litre,...
             'Ca+2 + H2O = CaOH+ + H+',  10^-12.78};

chem = ChemicalModel(elements, species, reactions);

chem.printChemicalSystem;

n = 100;

Na  = 1e-2*ones(n,1);
Cl  = 1e-2*ones(n,1);
H2O = ones(n,1);
Ca  = 1e-3*ones(n,1);
H   = logspace(-3, -11,n)';

inputs = [Na, Cl, Ca, H, H2O]*mol/litre;

tic
[state, report, chem] = chem.initState(inputs, 'chargeBalance', 'Na');
toc

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

addpath /Users/cmcneece/GoogleDrive/phreeqc/myfiles/comparisons/

DBname = 'phreeqc.dat';
filename = 'aqueousSpeciation';
progpath = '/Users/cmcneece/GoogleDrive/phreeqc/';

PHpath = ['/Users/cmcneece/GoogleDrive/phreeqc/myfiles/comparisons/'];
DBpath = '/Users/cmcneece/GoogleDrive/phreeqc/database/';
shellname =[PHpath filename];
    
    
    
eval(['! sh ' shellname '.sh ' shellname]);
eval(['! /Users/cmcneece/GoogleDrive/phreeqc/bin/phreeqc ' shellname '.txt ' shellname '.log ' DBpath DBname ]);

D = importdata([shellname '.sel']);

data = D.data;
pH = data(:,1);

ind = 2:size(data,2);
ind(3) =[];

plot(pH, data(:,ind),'--k', 'linewidth', 2)
xlim([3 11])
