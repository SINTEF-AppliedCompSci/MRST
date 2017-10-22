clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on


%% generate chemical system 

% define elements names
elements = {'C', 'O', 'H','Ca*'};

% define species names
species = {'CO2', 'HCO3-', 'CO3-2', 'H2CO3','Ca+2',...
           'H2O*','H+*','OH-',...
           'CaCO3(s)*', 'CO2(g)'};
        

% list chemical reactions         
reactions ={'CO2(g)  = CO2',              10^-1.47*(mol/litre)./atm,...
            'H2O = H+ + OH-',             10^-14*mol/litre,...
            'CO2 + H2O = H2CO3',          10^0./(mol/litre),...
            'H2CO3 = H+ + HCO3-',         10^-3.7*mol/litre,...
            'HCO3- = H+ + CO3-2',         10^-10.33*mol/litre,...
            'CaCO3(s) =  Ca+2 + CO3-2',   10^-7.*(mol/litre)^2};       

        
combinations = {'CT', 'CO2 + HCO3- + H2CO3 + CO3-2'};

% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions, 'combinations', combinations);

chem.plotIter = false;

% print the chemical system
chem.printChemicalSystem;

%% rock properties
n = 100;


%% solve the chemical system given inputs

H2O = ones(n,1)*mol/litre;
H = 1e-7*ones(n,1)*mol/litre;
Ca = 1e-1.*ones(n,1).*mol/litre;
CO2 = logspace(-2, 2, n)';


userInput = [Ca H2O H CO2];

tic
[state, report, model] = chem.initState(userInput);
toc;

[state, chem] = computeActivities(chem, state);
state = changeUnits(state, {'activities','masterComponents', 'components'}, mol/litre);

pH = -log10(chem.getProps(state, 'aH+'));
pH = log10(CO2);

 
CO2 = chem.getProps(state, 'CO2');
CO3 = chem.getProps(state, 'CO3-2');
HCO3 = chem.getProps(state, 'HCO3-');
H2CO3 = chem.getProps(state, 'H2CO3');


figure;hold on;
plot(pH, log10(state.components),'linewidth',2);
legend(chem.componentNames)


figure;hold on;
plot(pH, log10([state.saturationIndicies state.partialPressures]),'linewidth',2);
legend([chem.solidNames,chem.gasNames])


figure;hold on;
plot(pH, state.combinationComponents,'linewidth',2);
legend([chem.combinationNames])