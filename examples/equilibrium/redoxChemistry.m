clear; 
close all;

mrstModule add ad-core geochemistry
mrstVerbose on

elements = {'O', 'H', 'N*', 'e*','Na*'};

species = {'H+*', 'OH-', 'H2O*',...
            'NH4+', 'NO3-', 'e-','N2', 'NO2-','Na+'};

reactions ={'H2O  = H+  + OH-',                       10^-14*mol/litre,... 
            'NO3- + 10*H+ + 8*e- = NH4+ + 3*H2O',     10^119./(mol/litre)^15,...
            'NO3- + 2*H+ + 2*e- = NO2- + H2O',        10^28./(mol/litre)^3,...
            '2*NO3- + 12*H+ + 10*e- = N2 + 6*H2O',    10^(2*103)./(mol/litre)^17};
        
chem = ChemicalModel(elements, species, reactions);

chem.printChemicalSystem;

%%
n = 100;

N = 1e-3*ones(n,1);
e = logspace(-15, 5, n)';
H = 1e-7*ones(n,1);
H2O = ones(n,1);
Na = 1e-2*ones(n,1);

[state, report] = chem.initState([N e Na H H2O]*mol/litre);

%%
[state, chem] = chem.computeActivities(state);
[state, chem] = chem.computeChargeBalance(state);

state = changeUnits(state, {'species', 'activities'}, mol/litre);

pe = -log10(chem.getProp(state, 'ae-'));

figure;
plot(pe, state.species, 'linewidth', 2);
set(gca, 'yscale', 'log');
ylabel('concentration [mol/litre]');
xlabel('pe');
legend(chem.speciesNames);

