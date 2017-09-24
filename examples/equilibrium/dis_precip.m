clear;
close all;

%% masoud example
current = pwd;
cd('/Users/cmcneece/Downloads/redisspreciptransportupdate/')
run example.m
cd(current)


%% load adi and geochemistry module
run /Users/cmcneece/Documents/MATLAB/mrst-core/startup.m
mrstModule add geochemistry ad-core 
mrstVerbose on


%% generate chemical system 
k_B = 1; k_C = 0.67; omeg_B = 3; omeg_C = 2; phi_P = .4;

% define elements names
elements = {'Ba*','Ca*','SO4*'};

% define species names
species = {'Ba+2','Ca+2','SO4-2',...
            'BaSO4(s)','CaSO4(s)'};
        

% list chemical reactions         
reactions ={'CaSO4(s) = Ca+2 + SO4-2',       (k_C*(mol/litre)^2),...
            'BaSO4(s) = Ba+2 + SO4-2',       (k_B*(mol/litre)^2)};       

% list solid densities
solidDensities = {'CaSO4(s)', omeg_C*(mol/litre), 'BaSO4(s)',  omeg_B*(mol/litre)};

% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions);

chem.plotIter = false;

% print the chemical system
chem.printChemicalSystem;

%% rock properties
n = 100;
rock.perm = 1*darcy*ones(n, 1);
rock.poro = 0.4.*ones(n, 1);

%% solve the chemical system given inputs


Ba = linspace(0.2, 0.8, n)';
Ca = linspace(0.6, 0.2, n)';

SO4 = Ba+Ca;

userInput = [Ba Ca SO4]*mol/litre;

tic
[state, ~, ~] = chem.initState(userInput, 'solid', solidDensities, 'rock', rock);
toc;


figure(1);hold on;
plot(log10(state.components*litre/mol),'--','linewidth',2);
legend(chem.componentNames)


poros = [state.fluidVolumeFraction state.solidVolumeFractions];

figure(2);hold on;
plot(Ba.*100./SO4, fliplr(poros),'--', 'linewidth',2);
legend(fliplr(['fluid void fraction' chem.solidNames]))

