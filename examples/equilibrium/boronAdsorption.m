clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose off
%% generate chemical system 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Goldberg et al. Soil Sci. Soc. Am. J. 64:1356?1363 (2000) %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define master component names
elements = {'O', 'H', 'Na*', 'Cl*', 'B*'};

species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaCl', 'H2O*',...
             'H3BO3', 'H4BO4-',...
             '>SO-' '>SOH', '>SOH2+','>SOH3BO3-'};

reactions ={'H2O  <-> H+  + OH- ',              10^-14*mol/litre, ...
            'NaCl <-> Na+ + Cl-',               10^1*mol/litre,...
            'H3BO3 + H2O <-> H4BO4- + H+',      10^-9.2,...
            '>SOH + H+ <-> >SOH2+',             10^8.31/(mol/litre),...
            '>SOH <-> >SO- + H+',               10^-11.82*mol/litre,...
            '>SOH + H3BO3 <-> >SOH3BO3- + H+'    10^-8.2};
        
geometry = [2.31*site/(nano*meter)^2 59*meter^2/gram 200*gram/litre];
sInfo = {geometry, 'ccm', 1.06,     '>SO-',     -1,...
                                                       '>SOH',    0,...
                                                       '>SOH2+',   1,...
                                                       '>SOH3BO3-'  -1};
surfaces ={ '>SO', sInfo};

% instantiate chemical model
chem = ChemicalModel(elements, species, reactions, 'surf', surfaces);

% print the chemical system
chem.printChemicalSystem;

%% solve the chemical system given inputs
n = 200;

Na = 1e-1.*ones(n,1);
Cl = 1e-1.*ones(n,1);
B = 4.63e-4.*ones(n,1);
H2O = ones(n,1);
H = logspace(-7,-11, n)';

userInput = [Na Cl B H H2O]*mol/litre;

tic
state = [];
[state, report, model] = chem.initState(userInput);
toc;

[state, chem] = chem.computeActivities(state);
[state, chem] = chem.computeChargeBalance(state);


%% data from Goldberg

GB =   [7.145895	0.2151590
        7.259432	0.2151301
        7.277191	0.2502701
        7.399525	0.2578790
        7.662131	0.3296291
        7.977074	0.3937258
        8.029960    0.4517771
        8.380156	0.5540653
        8.669038	0.6349768
        9.097366	0.6807084
        9.176045	0.6898564
        9.604094	0.7019716
        9.813089	0.6285734
        10.29323	0.6040028
        10.65009	0.4572219
        10.81594	0.4464835];

%% take out relevant values from state

state = changeUnits(state, mol/litre );

SO      = getProp(chem, state, '>SO-');
SOH     = getProp(chem, state, '>SOH');
SOH2     = getProp(chem, state, '>SOH2+');

SOH3BO3     = getProp(chem, state, '>SOH3BO3-');

H      	= getProp(chem, state, 'H+');
OH     	= getProp(chem, state, 'OH-');
H3BO3	= getProp(chem, state, 'H3BO3');
H4BO4	= getProp(chem, state, 'H4BO4-');

pH = -log10(getProp(chem, state, 'aH+'));

%% phreeqc

chem = 'boron';

%% before stuff
addpath('/Users/cmcneece/GoogleDrive/phreeqc/');

% %% call and run phreeqc
% 
% current = pwd;
% 
% DBname = 'phreeqc.dat';
% filename = [chem '_sorption'];
% progpath = '/Users/cmcneece/GoogleDrive/phreeqc/';
% 
% PHpath = ['/Users/cmcneece/GoogleDrive/phreeqc/myfiles/' chem '/'];
% DBpath = '/Users/cmcneece/GoogleDrive/phreeqc/database/';
% shellname =[PHpath filename];
% 
% 
% 
% cd(progpath);
% eval(['! sh ' shellname '.sh ' shellname]);
% eval(['! ./bin/phreeqc ' shellname '.txt ' shellname '.log ' DBpath DBname ]);
% cd(current)
% 
% D = importdata([shellname '.sel']);
% 
% 
% 
% p.pH = D.data(:,1);
% p.S_B = D.data(:,5);
% p.H3BO3 = D.data(:,3);
% p.H4BO4 = D.data(:,4);
% p.H = D.data(:,6);
% p.OH = D.data(:,7);
% p.e = D.data(:,2); 

%% plot it
figure; hold on; box on;
plot(pH, state.chargebalance,'-k')
xlabel('pH')
ylabel('charge [% of total ion concentration]');

% aqueous components
figure; hold on; box on;
plot(pH, H3BO3, '-');
plot(pH, H4BO4, '-');
% plot(p.pH, p.H3BO3, '--r');
% plot(p.pH, p.H4BO4, '--r');
legend('H_3BO_3(aq)','4_3BO_4^-(aq)');

xlabel('pH');
ylabel('aqueous concentration [mol/L]');

% surface components
figure; hold on; box on;
plot(pH, SOH3BO3/200*1e6, '-k')
% plot(p.pH, p.S_B/200*1e6, '--r')
plot(GB(:,1), GB(:,2), 'ok');
legend('mrst', 'data')
   
xlabel('pH')
ylabel('adsorbed B [\mu mol/g]');

%%

