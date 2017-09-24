clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on

fromLoad = true;
%% generate chemical system 

% define master component names
elements = {'O', 'H', 'Ca*', 'C','Na*','Cl*'};

species = {'H+*', 'OH-','H2O*',...
            'Ca+2', 'CaCO3','CaHCO3+'...
            'CO3-2','HCO3-','CO2*',...
            'Na+', 'Cl-','NaCl'};

reactions ={'H2O  = H+  + OH- ',          10^-14*mol/litre,...
            'CO3-2 + 2*H+ = CO2 + H2O',   10^16.681/(mol/litre),...
            'HCO3- = H+ + CO3-2',      	10^-10.329*mol/litre,...
            'CaCO3 = Ca+2 + CO3-2',       10^3.224*mol/litre,...
            'Ca+2 + CO3-2 + H+ = CaHCO3+' 10^11.435/(mol/litre)^2,...
            'NaCl = Na+ + Cl-',           10^1*mol/litre};
                                                        
% instantiate chemical model
chem = ChemicalModel(elements, species, reactions);

% print the chemical system
chem.printChemicalSystem;

%% solve the chemical system given inputs
n = 100;

Ca = 1e-3.*ones(n,1);
Na = 1e-1.*ones(n,1);
Cl = 1e-1.*ones(n,1);
H = logspace(-1,-10, n)';
H2O = ones(n,1);
CO2 = 10^-1.468*10^-3.5.*ones(n,1);

userInput = [Ca Na Cl H H2O CO2]*mol/litre;

tic
[state, report, model] = chem.initState(userInput, 'charge', 'Cl');
toc;

[state, chem] = chem.computeActivities(state);
[state, chem] = chem.computeChargeBalance(state);
plot(log10(H),state.chargeBalance*litre/mol)


%% take out relevant values from state

state = changeUnits(state, {'masterComponents', 'components', 'activities'}, mol/litre );

pH = -log10(getProp(chem, state, 'aH+'));

%% phreeqc

folderName = 'mrstExamples';
filename ='carbonateSpeciation';

% call and run phreeqc

current = pwd;

DBname = 'phreeqc_simp.dat';
progpath = '/Users/cmcneece/GoogleDrive/phreeqc/';

PHpath = ['/Users/cmcneece/GoogleDrive/phreeqc/myfiles/' folderName '/'];
DBpath = '/Users/cmcneece/GoogleDrive/phreeqc/database/';
shellname =[PHpath filename];


if ~fromLoad
    cd(progpath);
    eval(['! sh ' shellname '.sh ' shellname]);
    eval(['! ./bin/phreeqc ' shellname '.txt ' shellname '.log ' DBpath DBname ]);
    cd(current)
end
D = importdata([shellname '.sel']);



p.pH = D.data(:,1);
p.H = D.data(:,4);
p.e = D.data(:,2); 
p.C = D.data(:,3);
p.Ca = D.data(:,8);
p.CaCO3 = D.data(:,9);

p.CO3 = D.data(:,5);
p.HCO3 = D.data(:,6);
p.CO2 = D.data(:,7);

 
 names = {'H+', 'C', 'CO2', 'HCO3-', 'CO3-2'};

figure; hold on; box on;
for i = 1 : numel(names)
    v = getProps(chem, state, names{i});
    plot(pH, v, '-')
end
plot(p.pH, p.H,'--k');
plot(p.pH, p.C,'--k');
plot(p.pH, p.CO3,'--k');
plot(p.pH, p.HCO3,'--k');
plot(p.pH, p.CO2,'--k');
xlabel('pH')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
legend(names)
ylim([1e-10 1e0]);
xlim([1 10]);
