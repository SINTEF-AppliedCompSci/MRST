clear;
close all;

%% load adi and geochemistry module
mrstModule add geochemistry ad-core 
mrstVerbose on


%% generate chemical system 
k_B = 1; k_C = 0.67;

% define elements names
elements = {'Ba*','Ca*','SO4'};

% define species names
species = {'Ba+2','Ca+2','SO4-2',...
            'BaSO4(s)*','CaSO4(s)'};
        

% list chemical reactions         
reactions ={'CaSO4(s) = Ca+2 + SO4-2',       (k_C*(mol/litre)^2),...
            'BaSO4(s) = Ba+2 + SO4-2',       (k_B*(mol/litre)^2)};       


% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions);

chem.plotIter = false;

% print the chemical system
chem.printChemicalSystem;

%% rock properties
n = 100;

%% solve the chemical system given inputs


Ba = logspace(-1, 1, n)'*mol/litre;
Ca = logspace(1, -1, n)'*mol/litre;

BaSO4 = 10*ones(n,1);

userInput = [Ba Ca BaSO4];

tic
[state, ~, ~] = chem.initState(userInput);
toc;


figure(1);hold on;
plot(log10(state.components*litre/mol),'--','linewidth',2);
legend(chem.componentNames)


figure(2);hold on;
plot(log10(state.saturationIndicies),'--','linewidth',2);
legend(chem.solidNames)
