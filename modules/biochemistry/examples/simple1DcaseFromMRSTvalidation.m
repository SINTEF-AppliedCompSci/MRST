%% Validation of MRST against two other simulators
% This example is a two-phase compositional problem in which CO2 is
% injected into a mixture of CO2, Methane, and Decane. The problem consists
% of 1000 cells and is one-dimensional for ease of visualization. The
% problem is divided into a large number of time steps to ensure that the
% different simulators take approximately the same timesteps.
%
% The problem is challenging in terms of fluid physics because the pressure
% is relatively low, which makes the phase behavior highly pressure
% dependent and all components exist in both phases. Since the wells are
% set to bottom-hole pressure controls, the fluid volume injected depends
% on correctly calculating the mobility and densities in the medium.
%
% MRST uses the Peng-Robinson equation of state by default and the Lohrenz,
% Bray, and Clark (LBC) correlation to determine viscosities for both
% phases.
%
% This example is discussed in Section 8.5.1 in the second MRST book:
% Advanced Modelling with the MATLAB Reservoir Simulation Toolbox (MRST),
% Cambridge University Press, 2021.
mrstModule add compositional deckformat ad-core ad-props

%% Set up model
% MRST includes both natural variables and overall composition. This toggle
% can switch between the modes.
if ~exist('useNatural', 'var')
    useNatural = true;
end
useNatural = false;
[state0, model, schedule, ref] = setupSimpleCompositionalExample(useNatural);
if useNatural
    name = 'Natural';
else
    name = 'Overall'; %#ok<UNRCH>
end
 eosname='pr';% 'pr';
%%
% Define synthetic saturation ranges
SW = linspace(0, 1, 100);  % Water saturation
SO = linspace(0, 1, 100);  % Oil saturation
SG = linspace(0, 1, 100);  % Gas saturation

% Evaluate original relative permeability functions
krW_original = arrayfun(model.fluid.krW, SW);  % Water rel perm
krOW_original = arrayfun(model.fluid.krOW, SO);  % Oil rel perm (water-oil)
krG_original = arrayfun(model.fluid.krG, SG);  % Gas rel perm
krOG_original = arrayfun(model.fluid.krOG, SO);  % Oil rel perm (gas-oil)

% Define residual saturations
swc = 0.2;  % Residual water saturation
sgr = 0.05;  % Residual gas saturation

% Update water relative permeability
SW_shifted = max(SW - swc, 0);  % Shift water saturation
krW_new = interp1(SW, krW_original, SW_shifted, 'linear', 0);  % Interpolate and set to 0 for sw <= swc
model.fluid.krW = @(sw) interp1(SW, krW_new, value(max(sw - swc, 0)));

% Update gas relative permeability
SG_shifted = max(SG - sgr, 0);  % Shift gas saturation
krG_new = interp1(SG, krG_original, SG_shifted, 'linear', 0);  % Interpolate and set to 0 for sg <= sgr
model.fluid.krG = @(sg) interp1(SG, krG_new, value(max(sg - sgr, 0)));

% Define residual oil saturation
sor = 0.1; % Example residual oil saturation

% Shift oil saturation for relative permeability
SO_shifted = max(SO - sor, 0); % Shift oil saturation
krOW_new = interp1(SO, krOW_original, SO_shifted, 'linear', 0); % Interpolate
krOG_new = interp1(SO, krOG_original, SO_shifted, 'linear', 0); % Interpolate

% Update oil relative permeability functions
model.fluid.krOW = @(so) interp1(SO, krOW_new, value(max(so - sor, 0)));
model.fluid.krOG = @(so) interp1(SO, krOG_new, value(max(so - sor, 0)));
model.fluid.krPts.w = [swc, 1];  % Water endpoints (residual and maximum)
model.fluid.krPts.g = [sgr, 1];  % Gas endpoints (residual and maximum)
% Define relative permeability endpoints
model.fluid.krPts.ow = [0, 1 - swc]; % Oil endpoints in water-oil system
model.fluid.krPts.og = [0, 1 - sgr]; % Oil endpoints in gas-oil system
%%% try dymanic 
% Define bacterial concentration parameters
phi0 = model.rock.poro; % Initial porosity
nc = 4e5; % Critical bacterial concentration (adjust as needed)
dp = 1e-5; % Pore diameter (adjust as needed)

% Initialize bacterial concentration in the state
state0.nbact = zeros(model.G.cells.num, 1); % Initial bacterial concentration

% Define porosity as a function of bacterial concentration
model.rock.poro = @(p, nbact) phi0 ./ (1 + 0.*(nbact / nc).^2);

pvMult_nbact = @(nbact) 1./ (1 + (0.*nbact / nc).^2);

% Replace the original pvMultR function
model.fluid.pvMultR = @(p, nbact) pvMult_nbact(nbact);

% Define permeability as a function of porosity

nbactMin = 1e3;
nbactMax = 4e5;
permMult = 1e-5;
tau = @(nbact) (min(max(nbact, nbactMin), nbactMax) - nbactMin)./(nbactMax - nbactMin);
perm = @(p,nbact) model.rock.perm(:,1).*((1-0.*tau(nbact)) + 0.*tau(nbact).*permMult);

model.rock.perm = perm;%@(p, nbact) (dp^2 / 180) * (model.rock.poro(p, nbact).^3) ./ (1 - model.rock.poro(p, nbact)).^2;

% Update the fluid object to use the new porosity and permeability
% model.fluid.poro = @(nbact) model.rock.poro(nbact);
% model.fluid.perm = @(nbact) model.rock.perm(nbact);

schedule.step.val = schedule.step.val/5;
 % Define compositional fluid model (with CoolProp library support)
 compFluid = TableCompositionalMixture({'Water', 'Hydrogen', 'CarbonDioxide', 'Methane'}, ...
     {'H2O', 'H2', 'CO2', 'C1'});
 model.EOSModel = SoreideWhitsonEquationOfStateModel(model.G, model.EOSModel.CompositionalMixture,eosname);
 diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true);
 mex_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'useMex', true, 'rowMajor', true);
 %includeWater=true
 arg = {model.G, model.rock, model.fluid, compFluid,true,diagonal_backend,...
     'water', false, 'oil', true, 'gas', true,'bacteriamodel', true,...
     'bDiffusionEffect', false,'moleculardiffusion',false,...
     'liquidPhase', 'O', 'vaporPhase', 'G'};
% model.setupOperators();
 model = BiochemistryModel(arg{:});
 model.outputFluxes = false;
 model.EOSModel.msalt=0;
 % Initial conditions
 z0 = [0.95, 1.0e-8, 0.1, 0.05-1.0e-8];
 schedule.control.W(1).components = [0, 0.95, 0.05, 0];
 schedule.control.W(2).components = [0, 0.95, 0.05, 0];
 T = 50 + 273.15;
 p = 75*barsa;
 nbact0 = 1.0e8;
 state0 = initCompositionalStateBacteria(model, p, T, [1, 0], z0, nbact0, model.EOSModel);
 model.verbose = false;
 problem = packSimulationProblem(state0, model, schedule, 'simple_comp_SW_bact_noclogging', 'name', name);
%% Simulate the schedule
% Note that as the problem has 500 control steps, this may take some time
% (upwards of 4 minutes).
simulatePackedProblem(problem, 'restartStep',1);
[ws, states, rep] = getPackedSimulatorOutput(problem);



%%
%% Compare with and without bectrial effects
problemNoBact = problem;
problemNoBact.BaseName = "simple_comp_SW_nobact_noclogging";
problemNoBact.OutputHandlers.states.dataDirectory= "/home/elyes/Documents/Projects/MRST/core/output/simple_comp_SW_nobact_noclogging";
problemNoBact.OutputHandlers.wellSols.dataDirectory= "/home/elyes/Documents/Projects/MRST/core/output/simple_comp_SW_nobact_noclogging";
problemNoBact.OutputHandlers.reports.dataDirectory= "/home/elyes/Documents/Projects/MRST/core/output/simple_comp_SW_nobact_noclogging";
[wsNoBact,statesNoBact] = getPackedSimulatorOutput(problemNoBact);
namecp = model.EOSModel.getComponentNames();
indH2=find(strcmp(namecp,'H2'));
indCO2= find(strcmp(namecp,'CO2'));
indCH4= find(strcmp(namecp,'C1'));
nT = numel(states);
% Initialize arrays to store total H2 mass
totalH2_bact = zeros(numel(states), 1);
totalH2_noBact = zeros(numel(statesNoBact), 1);

for i = 1:numel(states)
    % With bacterial effects
    totalH2_bact(i) = sum(states{i}.FlowProps.ComponentTotalMass{indH2});

    % Without bacterial effects
    totalH2_noBact(i) = sum(statesNoBact{i}.FlowProps.ComponentTotalMass{indH2});

    % With bacterial effects
    totalCO2_bact(i) = sum(states{i}.FlowProps.ComponentTotalMass{indCO2});

    % Without bacterial effects
    totalCO2_noBact(i) = sum(statesNoBact{i}.FlowProps.ComponentTotalMass{indCO2});
    % With bacterial effects
    totalCH4_bact(i) = sum(states{i}.FlowProps.ComponentTotalMass{indCH4});

    % Without bacterial effects
    totalCH4_noBact(i) = sum(statesNoBact{i}.FlowProps.ComponentTotalMass{indCH4});
end

%% Calculate percentage of H2 loss
H2_loss_percentage = ((totalH2_noBact - totalH2_bact) ./ totalH2_noBact) * 100;
%% Calculate percentage of CO2 loss
CO2_loss_percentage = ((totalCO2_noBact - totalCO2_bact) ./ totalCO2_noBact) * 100;
%% Calculate percentage of CH4 production
CH4_loss_percentage = ((totalCH4_bact - totalCH4_noBact) ./ totalCH4_noBact) * 100;

%% Display final H2 loss
fprintf('Total H2 loss due to bacterial effects: %.2f%%\n', H2_loss_percentage(end));
%% Comparison plots with existing simulators
% The same problem was defined into two other simulators: Eclipse 300
% (which is a commercial simulator) and AD-GPRS (Stanford's research
% simulator). We load in precomputed states from these simulators and
% compare the results.
%
% Note that the shock speed is sensitive to the different tolerances in the
% simulators, which have not been adjusted from the default in either
% simulator. We observe good agreement between all three simulators, with
% the minor differences can likely be accounted for by harmonizing the
% tolerances and sub-timestepping strategy for the different simulators.
ncomp = model.EOSModel.getNumberOfComponents();

lf = get(0, 'DefaultFigurePosition');
h = figure('Position', lf + [0, 0, 350, 0]);
data = {states};%, ref.statesECL(2:end), ref.statesGPRS};
n = min(cellfun(@numel, data));
names = {'MRST', 'E300', 'AD-GPRS'};
markers = {'-', '--', '--'};
cnames = model.EOSModel.getComponentNames();

nd = numel(data);
l = cell(nd*ncomp, 1);
for i = 1:nd
    for j = 1:ncomp
        l{(i-1)*ncomp + j} = [names{i}, ' ', cnames{j}];
    end
end
lw = [.5, 2.5, 5];
colors = lines(ncomp+1);
figure(h)
for step = 1:n % 180 for plot in book
    cla; hold on
    for i = 1:numel(data)
        s = data{i}{step};
        comp = s.components;
        nbact = s.nbact;
        if iscell(comp)
            comp = [comp{:}];
        end
        for j = 1:ncomp
            plot(comp(:, j), markers{i}, 'linewidth', lw(i), 'color', colors(j, :));
        end
        plot(nbact/max(states{end}.nbact) , markers{i}, 'linewidth', lw(2), 'color', colors(ncomp+1, :));

    end
    legend(l, 'location', 'north', 'numcolumns', 3);
    ylim([0, 1]);
    ylabel('z')
    drawnow
end

%% Compare pressure and saturations
% We also plot a more detailed comparison between MRST and E300 to show
% that the prediction of phase behavior is accurate.

colors = lines(ncomp + 2);
for step =  1:n
    figure(h); clf; hold on
    for i = 1:2
        s = data{i}{step};
        if i == 1
            marker = '-';
            linewidth = 1;
        else
            marker = '--';
            linewidth = 2.5;
        end
        hs = plot(s.s(:, 2), marker, 'color', [0.913, 0.172, 0.047], 'linewidth', linewidth, 'color', colors(1, :));
        p = s.pressure./max(s.pressure);
        hp = plot(p, marker, 'linewidth', linewidth, 'color', colors(2, :));
        comp = s.components;
        if iscell(comp)
            comp = [comp{:}];
        end
        
        if i == 1
            handles = [hs; hp];
        end
        for j = 1:ncomp
            htmp = plot(comp(:, j), marker, 'linewidth', linewidth, 'color', colors(j + 2, :));
            if i == 1
                handles = [handles; htmp]; %#ok<AGROW>
            end
        end
        if i == 2
            legend(handles, 'sV', 'Normalized pressure', cnames{:}, 'location', 'northoutside', 'orientation', 'horizontal');
        end
    end
    ylim([0, 1]);
    drawnow
end

%% Set up interactive plotting
% Finally we set up interactive plots to make it easy to look at the
% results from the different simulators.

mrstModule add mrst-gui
for i = 1:nd
    figure;
    plotToolbar(model.G, data{i}, 'plot1d', true);
    title(names{i});
end

%% Plot displacement lines in ternary diagram
figure; hold on
plot([0, 0.5, 1, 0], [0, sqrt(3)/2, 0, 0], 'k')


mapx = @(x, y, z) (1/2)*(2*y + z)./(x + y+ z);
mapy = @(x, y, z) (sqrt(3)/2)*z./(x + y+ z);

colors = parula(numel(states));
for i = 1:20:numel(states)
    C = states{i}.components;
    plot(mapx(C(:, 1), C(:, 2), C(:, 3)), mapy(C(:, 1), C(:, 2), C(:, 3)), '-', 'color', colors(i, :))
end
axis off

text(0, 0, cnames{1}, 'verticalalignment', 'top', 'horizontalalignment', 'right')
text(1, 0, cnames{2}, 'verticalalignment', 'top', 'horizontalalignment', 'left')
text(0.5, sqrt(3)/2, cnames{3}, 'verticalalignment', 'bottom', 'horizontalalignment', 'center')


text(mapx(0.5, 0.5, 0), mapy(0.5, 0.5, 0), '0.5', 'verticalalignment', 'top', 'horizontalalignment', 'center')
text(mapx(0, 0.5, 0.5), mapy(0, 0.5, 0.5), '0.5', 'verticalalignment', 'bottom', 'horizontalalignment', 'left')
text(mapx(0.5, 0.0, 0.5), mapy(0.5, 0.0, 0.5), '0.5', 'verticalalignment', 'bottom', 'horizontalalignment', 'right')

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
