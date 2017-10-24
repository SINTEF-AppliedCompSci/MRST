clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on

%% generate chemical system 

% define elements names
elements = {'O', 'H', 'Na*', 'Cl*'};

% define chemical species
species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaOH', 'H2O*',...
             '>SO-', '>SOH', '>SOH2+', '>SONa', '>SOH2Cl'};

% list chemical reactions         
reactions ={'H2O  = H+  + OH- ',          10^-14*mol/litre, ...
            'NaOH = Na+ + OH-',           10^10*mol/litre,...
            '>SOH = >SO- + H+',           10^-7.5*mol/litre,...
            '>SOH + H+ = >SOH2+',         10^3/(mol/litre),...
            '>SO- + Na+ = >SONa',         10^-2/(mol/litre),...
            '>SOH2+ + Cl- = >SOH2Cl',     10^-2/(mol/litre)};
        
% define the surface
geometry = [1*site/(nano*meter)^2 1*meter^2/gram 1*gram/litre];
sioInfo = {geometry, 'tlm', [1 0.2],	'>SONa',   [-1 1], '>SOH2Cl',[1 -1]};
surfaces ={ '>SO', sioInfo };

% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions, 'surf', surfaces);

% print the chemical system
chem.printChemicalSystem;

%%
n = 100;

Na = 1e-2*ones(n,1);
Cl = 1e-2*ones(n,1);
H = logspace(-3, -11, n)';
H2O = ones(n,1);

state = chem.initState([Na Cl H H2O]*mol/litre, 'chargeBalance', 'Na');

%%

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
filename = 'tripleLayerTest';
progpath = '/Users/cmcneece/GoogleDrive/phreeqc/';

PHpath = ['/Users/cmcneece/GoogleDrive/phreeqc/myfiles/comparisons/'];
DBpath = '/Users/cmcneece/GoogleDrive/phreeqc/database/';
shellname =[PHpath filename];
    
    
    
eval(['! sh ' shellname '.sh ' shellname]);
eval(['! /Users/cmcneece/GoogleDrive/phreeqc/bin/phreeqc ' shellname '.txt ' shellname '.log ' DBpath DBname ]);

D = importdata([shellname '.sel']);

data = D.data;
pH = data(:,1);

plot(pH, data(:,2:end),'--k', 'linewidth', 2)
xlim([3 11])