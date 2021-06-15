%% Comparison of the Natural Variable and Overall Composition Formulations
% This script contains examples from Section 8.4.3 of the second MRST book:
% Advanced Modelling with the MATLAB Reservoir Simulation Toolbox (MRST),
% Cambridge University Press, 2021.
%
% To demonstrates the difference in practice between the two compositional
% formulations implemented in MRST, and also show to set up the
% corresponding solvers, we consider a simplified two-component brine-CO2
% system posed inside a single grid cell.
mrstModule add compositional ad-core ad-props
mrstVerbose on

%% Construct the two simulation models
% Grid and petrophysical data in a single cell
G    = computeGeometry(cartGrid(1));      % Single cell 1 m^3 grid
rock = makeRock(G, 0.5*darcy, 0.5);

% Start with fluid properties from a standard black-oil, water-gas model
% and add extra mixture properties to turn this into a brine-CO2 model.
f = initSimpleADIFluid('phases', 'wg', 'blackoil', false, 'rho', [1000, 700]);
mixture = TableCompositionalMixture({'Water', 'CarbonDioxide'},{'Water', 'CO2'});

% Construct models for both formulations. Same input arguments
arg = {G, rock, f, ...                              % Standard arguments
       mixture,...                                  % Compositional mixture
       'water', true, 'oil', false, 'gas', true,... % Water-Gas system
       'liquidPhase', 'W', 'vaporPhase', 'G'};      % Water=liquid, gas=vapor
overall = GenericOverallCompositionModel(arg{:});   % Overall mole fractions
natural = GenericNaturalVariablesModel(arg{:});     % Natural variables

% Validate both models to initialize the necessary state function groups
overall = overall.validateModel();
natural = natural.validateModel();

%% Sidetrack: show state function groups for the overall model
% As a minor sidetrack, we produce plots of state functions like the ones
% used in the discussion of the implementation in Section 8.4.5

ogroups = overall.getStateFunctionGroupings();

% Plot all the end nodes for the FlowProperty group
% A reworked version of this plot appears in Figure 8.10.
figure
endnodes = {'CapillaryPressure', 'ComponentPhaseDensity',...
            'ComponentTotalMass','ComponentMobility'};
[~,g]=plotStateFunctionGroupings(ogroups(1:2),'Stop',endnodes,'label','name');
printStateFunctionGroupingTikz(g);
% printStateFunctionGroupingTikz(g,'file','compFlowProperty2.tex')

%  Plot all the end nodes for the PVTProperty group
figure
endnodes = {'ComponentPhaseMassFractions', 'Fugacity', 'Viscosity','ShrinkageFactors'};
plotStateFunctionGroupings(ogroups(1:2),'Stop',endnodes,'label','label');
% printStateFunctionGroupingTikz(g,'file','compPVTProperty.tex')

% Plot the entire graph
figure, plotStateFunctionGroupings(ogroups,'Stop', ...
    {'FlowDiscretization.ComponentTotalFlux',...
    'ComponentTotalMass',...
    'FacilityFlowDiscretization.ComponentTotalFlux'},'label','label');


%% Set up BC + initial conditions/initial guess
% We start by the boundary conditions. Here, we use standard routines from
% MRST, but must remember to also impose component specification at all
% boundaries
p = 50*barsa; T = 273.15 + 30; s = []; z = [1, 0];     % p, T, z, s
bc = fluxside([], G, 'xmin', 1/day, 'sat', [0, 1]);    % Flux
bc = pside(bc, G, 'xmax', p, 'sat', [0, 1]);           % Standard bc
bc.components = repmat([0, 1], numel(bc.face), 1);     % Boundary z

% We first set the initial state from the prescribed parameters, validate
% the model (which will add data structures for any wells, etc), and then
% perform a flash to ensure that the corresponding initial guess exists at
% equilibrium conditions
state0 = initCompositionalState(overall, p, T, s, z);
state0 = overall.validateState(state0);
initGuess = state0;
initGuess = overall.computeFlash(initGuess);

%% Adjust solver settings and construct nonlinear solver
% To improve the nonlinear solver we impose limits on changes during the
% nonlinear iterations. With natural variables, we limit phase mole
% fractions and saturation, whereas for the overall composition we limit
% the overall mole fraction
    % change
if ~exist('maxChange', 'var')
    maxChange = 0.1;
end
if isfinite(maxChange)
    overall.dzMaxAbs = maxChange;
    natural.dzMaxAbs = maxChange;
    natural.dsMaxAbs = maxChange;
end

% Set up nonlinear solver with high report level to output intermediate
% states
nls = NonLinearSolver('reportLevel', 3, 'MaxIterations', 1000);


%% Compute a single time step
dt = 100*day;
[solNat, reportsNat] = nls.solveTimestep(state0, dt, natural, ...
                            'bc', bc, 'initialGuess', initGuess);
[solMole, reportsMole] = nls.solveTimestep(state0, dt, overall, ...
                            'bc', bc, 'initialGuess', initGuess);

%% Extract data from the output states
getStates = @(reports) cellfun(@(x) x.state, reports.StepReports{1}.NonlinearReport, 'UniformOutput', false);
% Get the outputs
getZ = @(states) [initGuess.components(1); cellfun(@(x) x.components(:, 1), states)];
getS = @(states) [initGuess.s(:, 1); cellfun(@(x) x.s(:, 1), states)];
getP = @(states) [initGuess.pressure(1); cellfun(@(x) x.pressure, states)] - p;
natStates = getStates(reportsNat);
moleStates = getStates(reportsMole);

zn = getZ(natStates);
sn = getS(natStates);
dpn = getP(natStates);

zm = getZ(moleStates);
sm = getS(moleStates);
dpm = getP(moleStates);

%% Plot the convergence histories for both models
% We make two different plots of the convergence: in (z_w,dp) space and in
% (S_l,dp) space. These plots are imposed on top of a contour map of S.
ns = 50;
ps = linspace(min(min(dpn), min(dpm)), 1.1*max(max(dpn), max(dpm)), ns);
zs = linspace(0.01, 1, ns);

[DP, Z] = meshgrid(ps, zs);
zz = reshape(Z, [], 1);

[L, ~, ~, Z_L, Z_V] = standaloneFlash(DP(:) + p, T, [zz, 1 - zz], overall.EOSModel);
Lg = reshape(L, ns, ns);
S  = reshape(L.*Z_L./(L.*Z_L + (1-L).*Z_V), ns, ns);

for i = 1:2
    if i == 1
        xn = zn; xm = zm; XX = Z; l = 'z_{water}';
    else
        xn = sn; xm = sm; XX = S; l = 'S_L';       
    end
    figure(i); clf; hold on
    contourf(XX, DP, S, 10);
    c = lines(4);
    c1 = c(2, :);
    c2 = c(3, :);
    if maxChange < 0.1
        style = '-';
    else
        style = '-o';
    end
    h1 = plot(xn, dpn, style, 'color', c1, 'linewidth', 1.2, ...
        'MarkerFaceColor', 'r', 'MarkerSize', 8); hold on
    h2 = plot(xm, dpm, style, 'color', c2, 'linewidth', 1.2, ...
        'MarkerFaceColor', 'y');
    h3 = plot(xm(end), dpm(end), 'x', 'color', [1, 1, 1]*0.7, ...
        'MarkerSize', 10, 'linewidth', 1.2);
    legend([h1; h2; h3], 'Natural variables', 'Overall composition', ...
        'Solution', 'location', 'southeast')
    ylabel('\Delta p [Pa]')
    xlabel(l)

    cmap = interp1([1; 100], [.1, .1, 1; 0.9, 0.9, 1], linspace(1, 100, 100));
    colormap(flipud(cmap))
end

%%
%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}