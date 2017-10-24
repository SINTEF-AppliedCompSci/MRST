clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on

%% generate chemical system 

% define elements names
elements = {'O', 'H', 'Na*', 'Cl*'};

% define chemical species
species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaCl', 'H2O*',...
             '>SO-', '>SOH', '>SOH2+', '>SONa', '>SOH2Cl'};

% list chemical reactions         
reactions ={'H2O  = H+  + OH- ',          10^-14*mol/litre, ...
            'NaCl = Na+ + Cl-',           10^1*mol/litre,...
            '>SOH = >SO- + H+',           10^-7.5*mol/litre,...
            '>SOH + H+ = >SOH2+',         10^2/(mol/litre),...
            '>SO- + Na+ = >SONa',         10^-1.9/(mol/litre),...
            '>SOH2+ + Cl- = >SOH2Cl',     10^-1.9/(mol/litre)};
        
% define the surface
geometry = [2*site/(nano*meter)^2 50e-3*meter^2/gram 5e3*gram/litre];
sioInfo = {geometry, 'tlm', [1 0.2]/meter^2,      '>SONa',   [-1 1], '>SOH2Cl',[1 -1]};
surfaces ={ '>SO', sioInfo };

% define any linear combinations

% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions, 'surf', surfaces);

% print the chemical system
chem.printChemicalSystem;

%% solve the chemical system given inputs
n = 500;

Na = 1e-2.*ones(n,1);
Cl = 1e-2*ones(n,1);
H2O = ones(n,1);
H = logspace(-4, -10, n)';

userInput = [Na Cl H H2O]*mol/litre;


[state, report, model] = chem.initState(userInput, 'chargeBalance', 'Na');

%% post processing
[state, chem] = chem.computeActivities(state);
[state, chem] = chem.computeChargeBalance(state);
[state, chem] = chem.computeSurfaceCharges(state);
[state, chem] = chem.computeSurfacePotentials(state);
[state, chem] = chem.computeAqueousConcentrations(state);
[state, chem] = chem.computeSurfaceConcentrations(state);

state = changeUnits(state, {'elements', 'species', 'activities'}, mol/litre);

pH = -log10(getProp(chem, state, 'aH+'));

%% plot it
figure; hold on; box on;
plot(pH, state.chargeBalance,'-k')
xlabel('pH')
ylabel('charge [% of total ion concentration]');

figure; hold on; box on;
plot(pH, state.surfaceConcentrations, 'linewidth',2)
xlabel('pH')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
legend(chem.surfaceConcentrationNames);

figure; hold on; box on;
plot(pH, state.surfacePotentials, 'linewidth',2)
xlabel('pH')
ylabel('potential [volts]');
legend(chem.surfacePotentialNames);

figure; hold on; box on;
plot(pH, state.surfaceCharges, 'linewidth',2)
xlabel('pH')
ylabel('charge density [C/m^2]');
legend(chem.surfaceChargeNames);