clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on

%% generate chemical system 

elements = {'O', 'H', 'Na*', 'Cl*', 'Ca*', 'C'};

species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaOH', 'H2O*',...
             'Ca+2', 'CO3-2', 'HCO3-', 'CO2',...
             'CaCO3(s)', 'CO2(g)*', 'NaHCO3', 'CaCO3', 'CaHCO3+','CaOH+', 'NaCO3-'};

reactions ={'H2O  = H+  + OH- ',            10^-14*mol/litre, ...
            'NaOH = Na+ + OH-',             10^10*mol/litre,...
            'CaCO3(s) = CO3-2 + Ca+2',      10^-8.48*(mol/litre)^2,...
            'CO3-2 + H+ = HCO3-',           10^10.329/(mol/litre),...
            'CO3-2 + 2*H+ = CO2 + H2O',     10^16.681/(mol/litre),...
            'CO2(g) = CO2',                 10^-1.468*(mol/litre)/atm,...
            'Na+ + HCO3- = NaHCO3',         10^-0.25/(mol/litre),...
            'Na+ + CO3-2 = NaCO3-',         10^1.27/(mol/litre),...
            'Ca+2 + CO3-2 + H+ = CaHCO3+',  10^11.435/(mol/litre)^2,...
            'Ca+2 + CO3-2 = CaCO3',         10^3.224/(mol/litre),...
            'Ca+2 + H2O = CaOH+ + H+',      10^-12.78};

chem = ChemicalModel(elements, species, reactions);


%%
n = 100;

Na = 1e-1*ones(n,1)*mol/litre;
Cl = 1e-1*ones(n,1)*mol/litre;
Ca = 1e-3*ones(n,1)*mol/litre;
H = 1e-7*ones(n,1)*mol/litre;
H2O = ones(n,1)*mol/litre;
CO2 = logspace(-3,1,n)'*atm;

state = chem.initState([Na, Cl, Ca, H, H2O, CO2]);

%%

state = changeUnits(state, {'elements','species'}, mol/litre );

[state, chem] = chem.computeActivities(state);

figure(1); hold on; box on;
plot(CO2/atm, state.species, 'linewidth', 2)

xlabel('CO_2(g) [atm]')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');
legend(chem.speciesNames)


figure(2); hold on; box on;
plot(CO2/atm, state.saturationIndicies, 'linewidth', 2)

xlabel('CO_2(g) [atm]')
ylabel('saturation index [mol/L]');
set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');
legend(chem.solidNames)

%%
addpath /Users/cmcneece/GoogleDrive/phreeqc/myfiles/comparisons/

DBname = 'phreeqc.dat';
filename = 'phasesTest';
progpath = '/Users/cmcneece/GoogleDrive/phreeqc/';

PHpath = ['/Users/cmcneece/GoogleDrive/phreeqc/myfiles/comparisons/'];
DBpath = '/Users/cmcneece/GoogleDrive/phreeqc/database/';
shellname =[PHpath filename];
    
    
    
eval(['! sh ' shellname '.sh ' shellname]);
eval(['! /Users/cmcneece/GoogleDrive/phreeqc/bin/phreeqc ' shellname '.txt ' shellname '.log ' DBpath DBname ]);

D = importdata([shellname '.sel']);

data = D.data;

figure(1);
plot(10.^data(:,end-1), data(:,2:end-2),'.k', 'linewidth', 2)

figure(2);
plot(10.^data(:,end-1), 10.^data(:,end),'.k', 'linewidth', 2)


