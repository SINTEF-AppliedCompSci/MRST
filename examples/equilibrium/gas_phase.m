clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on


%% generate chemical system 

% define elements names
elements = {'C', 'O', 'H','Ca*'};

% define species names
species = {'CO2', 'CO2(g)*', 'HCO3-', 'CO3-2', 'H2O*', 'CaCO3(s)','H+*', 'Ca+2', 'H2CO3'};
% species = {'CO2', 'CO2(g)*', 'HCO3-', 'CO3-2', 'H2O*','H+*', 'Ca+2', 'H2CO3'};
        

% list chemical reactions         
reactions ={'CO2(g)  = CO2',              10^-1.47*(mol/litre)./atm,...
            'CO2 + H2O = H2CO3',          10^0./(mol/litre),...
            'H2CO3 = H+ + HCO3-',         10^-3.7*mol/litre,...
            'HCO3- = H+ + CO3-2',         10^-10.33*mol/litre,...
            'CaCO3(s) =  Ca+2 + CO3-2',   10^-6.36.*mol/litre};       

% reactions ={'CO2(g)  = CO2',              10^-1.47*(mol/litre)./atm,...
%             'CO2 + H2O = H2CO3',          10^-1.9./(mol/litre),...
%             'H2CO3 = H+ + HCO3-',         10^-3.7*mol/litre,...
%             'HCO3- = H+ + CO3-2',         10^-10.33*mol/litre};  
        
combinations = {'CT', 'CO2 + HCO3- + CO3-2 + H2CO3'};

% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions); %, 'combinations', combinations);

chem.plotIter = false;

% print the chemical system
chem.printChemicalSystem;

%% rock properties
n = 100;
rock.perm = 1*darcy*ones(n, 1);
rock.poro = 0.4.*ones(n, 1);

%% solve the chemical system given inputs

CT = 1e-1.*ones(n,1);
H2O = ones(n,1)*mol/litre;
% H = 1e-7.*ones(n,1).*mol/litre;
H = logspace(-7, -10,n)'.*mol/litre;
Ca = 1e1.*ones(n,1).*mol/litre;
% CO2 = logspace(-5, 1,n)'.*atm;
CO2 = 10^-1.54.*ones(n,1);

% list solid densities
solidDensities = {'CaCO3(s)', 3*mol/litre};

% list partial pressures
partialPressure = {'CO2(g)', 1*atm};

userInput = [Ca CO2 H2O H];
% userInput = [C O H Ca]*mol/litre;

tic
[state, report, model] = chem.initState(userInput, 'solid', solidDensities);
% [state, report, model] = chem.initState(userInput, 'partial', partialPressure);
toc;
[state, chem] = computeActivities(chem, state);
state = changeUnits(state, {'activities','masterComponents', 'components'}, mol/litre);

pH = -log10(chem.getProps(state, 'aH+'));

CO2 = chem.getProps(state, 'CO2');
CO3 = chem.getProps(state, 'CO3-2');
HCO3 = chem.getProps(state, 'HCO3-');
H2CO3 = chem.getProps(state, 'H2CO3');

CT = CO2 + CO3 + HCO3 + H2CO3;


figure;hold on;
plot(pH, log10(state.components),'linewidth',2);
plot(pH, log10(CT),'linewidth',2);
legend([chem.CompNames 'DIC'])


figure;hold on;
plot(pH, [state.solidComponents state.gasComponents state.poro],'linewidth',2);
legend([chem.SolidNames,chem.GasNames 'porosity'])


