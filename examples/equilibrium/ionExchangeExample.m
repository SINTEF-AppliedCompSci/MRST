clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on

fromLoad = true;

%% generate chemical system 

% define elements names
elements = {'O', 'H', 'Ca*', 'Na*','Cl*'};

% define chemical species
species = {'H+*', 'OH-', 'Ca+2','Cl-', 'Na+', 'H2O*',...
             '>XNa', '>X2Ca', '>XH'};

% list chemical reactions         
reactions ={'H2O  = H+  + OH- ',              10^-14*mol/litre, ...
            '>XH + Na+ = >XNa + H+',       	10^-1,...
            '2*>XH + Ca+2 = >X2Ca + 2*H+',  	10^-1.2};
        
% define the surface
geometry = [2*site/(nano*meter)^2 50e-3*meter^2/(gram) 5e3*gram/litre];
xInfo = {geometry, 'ie'};

surfaces ={ '>X', xInfo };
                                                        
% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions, 'surf', surfaces);

chem.plotIter = false;

% print the chemical system
chem.printChemicalSystem;

%% solve the chemical system given inputs
n = 500;

Ca = logspace(-4,-1,n)';
Na = 1e-2*ones(n,1);
Cl = 3e-2.*ones(n,1);
H2O = ones(n,1);
H = 1e-4.*ones(n,1);

userInput = [Ca Na Cl H H2O]*mol/litre;

tic
[state, report, model] = chem.initState(userInput, 'chargeBalance', 'Cl');
toc;

[state, chem] = chem.computeActivities(state);
[state, chem] = chem.computeChargeBalance(state);


%% plot

state = changeUnits(state, {'masterComponents', 'components'}, mol/litre );

Ca = getProp(chem, state, 'Ca+2');

figure; hold on; box on;

plot(Ca, state.components, 'linewidth', 2)
set(gca, 'yscale','log')
set(gca, 'xscale','log')
xlabel('Ca / M')
ylabel('concentration / M');
xlim([min(Ca), max(Ca)]);
legend(chem.componentNames);
