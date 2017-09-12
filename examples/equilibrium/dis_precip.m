clear;
close all;

% load adi and geochemistry module
mrstModule add geochemistry ad-core 
mrstVerbose on


%% generate chemical system 

% define elements names
elements = {'Ba*','Ca*','SO4*'};

% define species names
species = {'Ba+2','Ca+2','SO4-2',...
            'BaSO4(s)','CaSO4(s)'};
        

% list chemical reactions         
reactions ={'CaSO4(s)  <-> Ca+2 + SO4-2 ',       0.67*mol/litre,...
            'BaSO4(s)  <-> Ba+2 + SO4-2',        1*mol/litre};       

% list solid densities
solidDensities = {'CaSO4(s)', 2*mol/litre, 'BaSO4(s)',  3*mol/litre};

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


Ba  =logspace(-1, 1,n)';
Ca  = logspace(1, -1,n)';

% Ba  =logspace(-2, 2,n)';
% Ca  = logspace(2, -2,n)';
SO4 = Ba+Ca;

userInput = [Ba Ca SO4]*mol/litre;
% userInput = [Ba Ca C H H2O]*mol/litre;

tic
[state, report, model] = chem.initState(userInput, 'solid', solidDensities);
toc;


figure;hold on;
plot(log10(state.components*litre/mol),'linewidth',2);
legend(chem.CompNames)


poros = [state.solidComponents state.poro];
for i = 1 : size(poros,2);
    poros(:,i) = poros(:,i).*rock.poro;
end

figure;hold on;
plot([rock.poro poros],'linewidth',2);
legend(['rock porosity', chem.SolidNames,'fluid porosity'])

