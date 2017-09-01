clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on


%% generate chemical system 

% define elements names
elements = {'C', 'O', 'H','Ca*'};

% define species names
species = {'CO2', 'CO2(g)', 'HCO3-', 'CO3-2', 'H2O*', 'CaCO3(s)','H+*', 'Ca+2', 'H2CO3'};
        

% list chemical reactions         
reactions ={'CO2(g)  <-> CO2',              10^-1.5*mol/litre./atm,...
            'CO2 + H2O <-> H2CO3',          10^-1.47./(mol/litre),...
            'H2CO3 <-> H+ + HCO3-',         10^-6.35*mol/litre,...
            'HCO3- <-> H+ + CO3-2',         10^-10.33*mol/litre,...
            'CaCO3(s) <->  Ca+2 + CO3-2',   10^-6.36.*mol/litre};       

combinations = {'CT*', 'CO2 + HCO3- + CO3-2 + H2CO3'};

% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions, 'combinations', combinations);

chem.plotIter = false;

% print the chemical system
chem.printChemicalSystem;

%% rock properties
n = 100;
rock.perm = 1*darcy*ones(n, 1);
rock.poro = 0.4.*ones(n, 1);

%% solve the chemical system given inputs

CT = 1e1.*ones(n,1);
H2O = ones(n,1);
H = logspace(-4, -10,n)';
Ca = 1e1.*ones(n,1);


% list solid densities
solidDensities = {'CaCO3(s)', 3*mol/litre};

% list partial pressures
partialPressure = {'CO2(g)', 1*atm};

userInput = [Ca H2O H CT]*mol/litre;
% userInput = [Ba Ca C H H2O]*mol/litre;

tic
[state, report, model] = chem.initState(userInput, 'solid', solidDensities,...
                                        'partial', partialPressure);
toc;


figure;hold on;
plot(log10(state.components*litre/mol),'linewidth',2);
legend(chem.CompNames)

figure;hold on;
plot([state.solidComponents state.poro],'linewidth',2);
legend([chem.SolidNames,'porosity'])


