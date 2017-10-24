clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on

%% generate chemical system 

elements = {'O', 'H', 'Na*', 'Cl*', 'Ca*', 'C'};

species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaOH', 'H2O*',...
             'Ca+2', 'CO3-2', 'HCO3-', 'CO2',...
             'CaCO3(s)', 'CO2(g)*'};

reactions ={'H2O  = H+  + OH- ',            10^-14*mol/litre, ...
            'NaOH = Na+ + OH-',             10^10*mol/litre,...
            'CaCO3(s) = CO3-2 + Ca+2',      10^-8.48*(mol/litre)^2,...
            'CO3-2 + H+ = HCO3-',           10^10.329/(mol/litre),...
            'CO3-2 + 2*H+ = CO2 + H2O',     10^16.681/(mol/litre),...
            'CO2 = CO2(g)',                 10^1.468*atm/(mol/litre)};

chem = ChemicalModel(elements, species, reactions);


%%
n = 100;

Na = 1e-2*ones(n,1)*mol/litre;
Cl = 1e-2*ones(n,1)*mol/litre;
H = 1e-7*ones(n,1)*mol/litre;
H2O = ones(n,1)*mol/litre;
Ca = 1e-3*ones(n,1)*mol/litre;
CO2 = logspace(-3, 1, n)'*atm;

state = chem.initState([Na, Cl, Ca, H, H2O, CO2]);

%%

state = changeUnits(state, {'elements','species','partialPressures'}, [mol/litre, mol/litre, atm] );

[state, chem] = chem.computeActivities(state);

figure; hold on; box on;
plot(CO2/atm, state.species, 'linewidth', 2)
set(gca, 'xscale','log');

xlabel('CO_2(g) [atm]')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
legend(chem.speciesNames)

figure; hold on; box on;
plot(CO2/atm, state.saturationIndicies, 'linewidth', 2)
set(gca, 'xscale','log');

xlabel('CO_2(g) [atm]')
ylabel('saturation indicies [-]');
set(gca, 'yscale', 'log');
legend(chem.solidNames)
