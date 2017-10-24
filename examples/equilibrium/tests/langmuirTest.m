%% create and solve chemical system
clear;
close all;

mrstModule add geochemistry ad-core
mrstVerbose on
%% generate chemical system object

elements = {'O', 'H', 'Na*', 'Cl*'};

species = {'H+*', 'OH-', 'H2O*', '>SiO-', '>SiOH', '>SiOH2+', 'Na+', 'Cl-','NaOH'};

reactions = {'H2O  = H+  + OH- ',           10^-14*mol/litre, ...
             '>SiOH = >SiO- + H+',          10^-7.5*mol/litre,...
             '>SiOH + H+ = >SiOH2+',         10^3/(mol/litre),...
             'NaOH = Na+ + OH-',            10^10*mol/litre};

surfInfo = {'>SiO', {[1*site/(nano*meter)^2 1*meter^2/gram 1*gram/litre], 'langmuir'}};

chem = ChemicalModel(elements, species, reactions, 'surfaces', surfInfo);


n = 100;

H2O = ones(n,1);
H   = logspace(-3, -11,n)';
Na  = 1e-2*ones(n,1);
Cl  = Na;

inputs = [Na, Cl, H, H2O]*mol/litre;

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
filename = 'langmuirTest';
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

plot(pH, data(:,2:end),'--k', 'linewidth', 2)
xlim([3 11])
