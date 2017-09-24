clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on

fromLoad = false;

%% generate chemical system 

% define elements names
elements = {'O', 'H', 'Na*', 'Cl*'};

% define chemical species
species = {'H+*', 'OH-', 'Na+', 'Cl-', 'NaCl', 'H2O*',...
             '>SO-', '>SOH', '>SOH2+', '>SONa', '>SOH2Cl'};

% list chemical reactions         
reactions ={'H2O  = H+  + OH- ',          10^-14*mol/litre, ...
            'NaCl = Na+ + Cl-',           10^1*mol/litre,...
            '>SOH = >SO- + H+',           10^-7.5*mol/litre,...
            '>SOH + H+ = >SOH2+',         10^2/(mol/litre),...
            '>SO- + Na+ = >SONa',         10^-1.9/(mol/litre),...
            '>SOH2+ + Cl- = >SOH2Cl',     10^-1.9/(mol/litre)};
        
% define the surface
geometry = [2*site/(nano*meter)^2 50e-3*meter^2/gram 5e3*gram/litre];
sioInfo = {geometry, 'tlm', [1 0.2]/meter^2,      '>SONa',   [-1 1], '>SOH2Cl',[1 -1]};
surfaces ={ '>SO', sioInfo };

% define any linear combinations

% instantiate the chemical model
chem = ChemicalModel(elements, species, reactions, 'surf', surfaces);

chem.plotIter = false;

% print the chemical system
chem.printChemicalSystem;

%% solve the chemical system given inputs
n = 500;

Na = 1e-1.*ones(n,1);
Cl = 1e-1*ones(n,1);
% B = 4.63e-4;
H2O = ones(n,1);
H = logspace(-4, -10, n)';

userInput = [Na Cl H H2O]*mol/litre;

tic
[state, report, model] = chem.initState(userInput, 'chargeBalance', 'Na');
toc;

[state, chem] = chem.computeActivities(state);
[state, chem] = chem.computeChargeBalance(state);
[state, chem] = chem.computeSurfaceCharges(state);
[state, chem] = chem.computeSurfacePotentials(state);


%% take out relevant values from state
state = changeUnits(state, {'masterComponents', 'components', 'activities'}, mol/litre);

SO      = getProp(chem, state, '>SO-');
SOH     = getProp(chem, state, '>SOH');
SOH2    = getProp(chem, state, '>SOH2+');

SOH2Cl     = getProp(chem, state, '>SOH2Cl');
SONa     = getProp(chem, state, '>SONa');

H      	= getProp(chem, state, 'H+');
OH     	= getProp(chem, state, 'OH-');

pH = -log10(getProp(chem, state, 'aH+'));

%% phreeqc

folderName = 'mrstExamples';
filename ='tripleLayerModelExample';

% call and run phreeqc

current = pwd;

DBname = 'phreeqc.dat';
progpath = '/Users/cmcneece/GoogleDrive/phreeqc/';

PHpath = ['/Users/cmcneece/GoogleDrive/phreeqc/myfiles/' folderName '/'];
DBpath = '/Users/cmcneece/GoogleDrive/phreeqc/database/';
shellname =[PHpath filename];


if ~fromLoad
    cd(progpath);
    eval(['! sh ' shellname '.sh ' shellname]);
    eval(['! ./bin/phreeqc ' shellname '.txt ' shellname '.log ' DBpath DBname  '&>/dev/null']);
    cd(current)
end

D = importdata([shellname '.sel']);

p.pH = D.data(:,1);
p.e = D.data(:,2); 
p.SO = D.data(:,3);
p.SOH = D.data(:,4);
p.SOH2 = D.data(:,5);
p.SONa = D.data(:,6);
p.SOH2Cl = D.data(:,7);




%% plot it
figure; hold on; box on;
plot(pH, state.chargeBalance,'-k')
plot(p.pH, p.e, '--r')
xlabel('pH')
ylabel('charge [% of total ion concentration]');
legend('mrst', 'phreeqc');

% surface components
figure; hold on; box on;
plot(pH, SO)
plot(pH, SOH)
plot(pH, SOH2)
plot(pH, SONa)
plot(pH, SOH2Cl)
plot(p.pH, p.SO, '--k')
plot(p.pH, p.SOH, '--k')
plot(p.pH, p.SOH2, '--k')
plot(p.pH, p.SONa, '--k')
plot(p.pH, p.SOH2Cl, '--k') 

xlim([4 10])
xlabel('pH')
ylabel('concentration [mol/L]');
set(gca, 'yscale', 'log');
legend('>SiO^-','>SiOH','>SiOH2^+','>SiONa','>SiOH_2Cl', 'location','southwest');
