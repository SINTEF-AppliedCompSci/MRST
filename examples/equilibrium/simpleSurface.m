clear;
close all;

% load adi and geochemistry module
mrstModule add ad-core geochemistry
mrstVerbose on

%% generate chemical system 

elements = {'O', 'H', 'Na*', 'Cl*'};

species = {'H+*', 'OH-', 'Na+', 'Cl-', 'H2O*', '>XH', '>XNa'};

reactions ={'H2O  = H+  + OH- ',          10^-14*mol/litre, ...
            '>XH + Na+ = >XNa + H+',    10^-8};
        
geometry = [1*site/(nano*meter)^2 50*meter^2/gram 1e3*gram/litre];
sioInfo = {geometry, 'ie',};
surfaces ={ '>X', sioInfo };

chem = ChemicalModel(elements, species, reactions, 'surf', surfaces);


% print the chemical system
chem.printChemicalSystem;

%% solve the chemical system given inputs
n = 50;l = 50;

H = logspace(-4, -10, n)';
Hrep = repmat(H,l,1);

Na = logspace(-4,-1,l);
Narep = repmat(Na,n,1);
Narep = Narep(:);

Cl = Narep;
H2O = ones(n*l,1);

userInput = [Narep Cl Hrep H2O]*mol/litre;

[state, report, model] = chem.initState(userInput);

state = changeUnits( state, {'species','elements'}, mol/litre);

[state, chem] = chem.computeActivities(state);
[state, chem] = chem.computeSurfaceConcentrations(state);
[state, chem] = chem.computeAqueousConcentrations(state);

[H_aq, Na_aq, H_surf, Na_surf] = chem.getProps(state, 'H(aq)', 'Na(aq)', 'H(surf)', 'Na(surf)');

H_surf = reshape(H_surf, n, l);
Na_surf = reshape(Na_surf, n, l);

figure;
subplot 121
contourf(H, Na, H_surf', 'linestyle','none');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('[H^+] / M');
ylabel('[Na^+] / M');
h = colorbar;
ylabel(h, '[XH] / M');
colormap parula

subplot 122
contourf(H, Na, Na_surf', 'linestyle','none');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('[H^+] / M');
ylabel('[Na^+] / M');
h = colorbar;
ylabel(h, '[XNa] / M');
colormap parula
