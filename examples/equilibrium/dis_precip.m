clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on


%% generate chemical system 

% define elements names
elements = {'Ba*','Ca*','SO4*'};

% define species names
species = {'Ba+2','Ca+2','SO4-2',...
            'BaSO4(s)','CaSO4(s)'};
        

% list chemical reactions         
reactions ={'CaSO4(s)  <-> Ca+2 + SO4-2 ',       1*mol/litre,...
            'BaSO4(s)  <-> Ba+2 + SO4-2',        0.67*mol/litre};       

% list solid densities
solidDensities = {'CaSO4(s)', 3*mol/litre, 'BaSO4(s)',  2*mol/litre};

% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions, 'solid', solidDensities);

chem.plotIter = false;

% print the chemical system
chem.printChemicalSystem;

%% rock properties
n = 10;
rock.perm = 1*darcy*ones(n, 1);
rock.poro = 0.4.*ones(n, 1);

%% solve the chemical system given inputs


Ba  =logspace(-1, 1,n)';
Ca  = logspace(-1, 1,n)';
SO4 = Ba + Ca;

userInput = [Ba Ca SO4]*mol/litre;
% userInput = [Ba Ca C H H2O]*mol/litre;

state.poro = rock.poro;
tic
[state, report, model] = chem.initState(userInput, 'state', state);
toc;

[state, chem] = chem.computeActivities(state);
[state, chem] = chem.computeChargeBalance(state);

figure;hold on;
plot(log10(CaCO3), log10(state.components*litre/mol),'linewidth',2);
legend(chem.CompNames)



%% call and run phreeqc
folderName = 'mrstExamples';
filename ='pure_phase_eq';

current = pwd;

DBname = 'dis_precip.dat';
progpath = '/Users/cmcneece/GoogleDrive/phreeqc/';

PHpath = ['/Users/cmcneece/GoogleDrive/phreeqc/myfiles/' folderName '/'];
DBpath = '/Users/cmcneece/GoogleDrive/phreeqc/database/';
shellname =[PHpath filename];



    cd(progpath);
%     eval(['! sh ' shellname '.sh ' shellname]);
    eval(['! ./bin/phreeqc ' shellname '.txt ' shellname '.log ' DBpath DBname ]);
    cd(current)

D = importdata([shellname '.sel']);


for i = 2 : 7
    plot((D.data(:,end-1)), log10(D.data(:,i)),'.k')
end

% for i = 3 : size(D.data,2)
%     plot(log10(D.data(:,2)), log10(D.data(:,i)),'.k')
% end

