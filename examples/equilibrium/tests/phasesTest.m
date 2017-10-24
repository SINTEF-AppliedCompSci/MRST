clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on

%% generate chemical system 

elements = {'O', 'H', 'Na*', 'Cl*', 'Ca', 'C'};

species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaOH', 'H2O*',...
             'Ca+2', 'CO3-2', 'HCO3-', 'CO2',...
             'CaCO3(s)*', 'CO2(g)*', 'NaHCO3', 'CaCO3', 'CaHCO3+','CaOH+', 'NaCO3-'};

reactions ={'H2O  = H+  + OH- ',            10^-14*mol/litre, ...
            'NaOH = Na+ + OH-',             10^10*mol/litre,...
            'CaCO3(s) = CO3-2 + Ca+2',      10^-8.48*(mol/litre)^2,...
            'CO3-2 + H+ = HCO3-',           10^10.329/(mol/litre),...
            'CO3-2 + 2*H+ = CO2 + H2O',     10^16.681/(mol/litre),...
            'CO2 = CO2(g)',                 10^1.468*atm/(mol/litre),...
            'Na+ + HCO3- = NaHCO3',         10^-0.25/(mol/litre),...
            'Na+ + CO3-2 = NaCO3-',         10^1.27/(mol/litre),...
            'Ca+2 + CO3-2 + H+ = CaHCO3+',  10^11.435/(mol/litre)^2,...
            'Ca+2 + CO3-2 = CaCO3',         10^3.224/(mol/litre),...
            'Ca+2 + H2O = CaOH+ + H+',      10^-12.78};

chem = ChemicalModel(elements, species, reactions);


%%
n = 100;

Na = 1e-2*ones(n,1)*mol/litre;
Cl = 1e-2*ones(n,1)*mol/litre;
H = 1e-7*ones(n,1)*mol/litre;
H2O = ones(n,1)*mol/litre;
CaCO3 = logspace(-1, 1, n)';
CO2 = 1e-1*ones(n,1)*atm;

state = chem.initState([Na, Cl, H, H2O, CaCO3, CO2]);

%%

state = changeUnits(state, {'elements','species'}, mol/litre );

[state, chem] = chem.computeActivities(state);

figure; hold on; box on;
plot(log10(CaCO3), state.species, 'linewidth', 2)

xlabel('CaCO3 log_{10}(saturation index)')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
legend(chem.speciesNames)

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

plot(data(:,end), data(:,2:end-2),'--k', 'linewidth', 2)
