clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on


%% generate chemical system 

% define elements names
elements = {'H','Ba','Ca','C*','O'};
% elements = {'H','Ba*','Ca*','C*','O'};

% define chemical species
species = {'H+*','OH-','H2O*' 'HCO3-',...
            'Ba+2','Ca+2','CO3-2',...
            'BaCO3','CaCO3',...
            'BaCO3(s)*','CaCO3(s)*'};
        
% species = {'H+*','OH-','H2O*' 'HCO3-',...
%             'Ba+2','Ca+2','CO3-2',...
%             'BaCO3','CaCO3'};
        
% list chemical reactions         
reactions ={'CO3-2 + H+ <-> HCO3-',         10^10.329/(mol/litre),...
            'H2O <-> H+ + OH-',              10^-14*mol/litre,...
            'CaCO3  <-> Ca+2 + CO3-2 ',       10.^-3.224*mol/litre,...
            'BaCO3  <-> Ba+2 + CO3-2',        10.^-2.71*mol/litre,...
            'CaCO3(s)  <-> Ca+2 + CO3-2 ',       10.^-8.48*mol/litre,...
            'BaCO3(s)  <-> Ba+2 + CO3-2',        10.^-8.562*mol/litre};
        
%  % list chemical reactions         
% reactions ={'CO3-2 + H+ <-> HCO3-',         10^10.329/(mol/litre),...
%             'H2O <-> H+ + OH-',              10^-14*mol/litre,...
%             'CaCO3  <-> Ca+2 + CO3-2 ',       10.^-3.224*mol/litre,...
%             'BaCO3  <-> Ba+2 + CO3-2',        10.^-2.71*mol/litre};       

% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions);

chem.plotIter = false;

% print the chemical system
chem.printChemicalSystem;

%% solve the chemical system given inputs
n = 100;
C = 1e-3.*ones(n,1);
H = 1e-7.*ones(n,1);
H2O = ones(n,1);
Ba = logspace(-1,0,n)';
Ca = 1e-3*ones(n,1);

BaCO3 = 1e-3*ones(n,1);
CaCO3 = logspace(-6,-3,n)';

userInput = [C H H2O BaCO3 CaCO3]*mol/litre;
% userInput = [Ba Ca C H H2O]*mol/litre;

state.porosity = 
tic
[state, report, model] = chem.initState(userInput, 'state',state);
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

