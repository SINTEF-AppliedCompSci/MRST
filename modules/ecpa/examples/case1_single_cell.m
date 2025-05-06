%% A two-component water-CO2 system inside a single grid cell
% This script contains the example case 1 for the paper
% "Compositional simulation for carbon storage in porous media using
%  an electrolyte association equation of state"
%
% To demonstrates the difference in practice between the two compositional
% formulations (Natural Variable and Overall Composition Formulations)
% implemented in MRST, and also show to set up the corresponding solvers,
% we consider a simplified two-component water-CO2 system posed inside
% a single grid cell.
mrstModule add compositional ad-core ad-props ecpa
mrstVerbose on

%% Construct the two simulation models
% Grid and petrophysical data in a single cell
G    = computeGeometry(cartGrid(1));      % Single grid cell with size 1 m^3
rock = makeRock(G, 0.5*darcy, 0.5);

% Start with fluid properties from a standard black-oil, water-gas model
% and add extra mixture properties to turn this into a brine-CO2 model.
f = initSimpleADIFluid('phases', 'wg', 'blackoil', false, 'rho', [1000, 700]);
T = 273.15 + 30;
mixture = ECPATableCompositionalMixture({'Water', 'CarbonDioxide'});
bic = eCPAreadBinaryInteraction(mixture,T(1));
ap = eCPAreadAssociationParameter(mixture, T(1));
mixture = setBinaryInteraction(mixture, bic);
mixture = setAssociationParameter(mixture, ap);

% Construct models for both formulations using the same input arguments
arg = {G, rock, f, ...                                  % Standard arguments
       mixture,...                                      % Compositional mixture
       'water', true, 'oil', false, 'gas', true,...     % Water-Gas system
       'liquidPhase', 'W', 'vaporPhase', 'G'};          % Water=liquid, gas=vapor
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', true, 'rowMajor', true);
overall = ECPAGenericOverallCompositionModel(arg{:}, 'AutoDiffBackend', diagonal_backend);   % Overall mole fractions
natural = ECPAGenericNaturalVariablesModel(arg{:}, 'AutoDiffBackend', diagonal_backend);     % Natural variables

% Validate both models to initialize the necessary state function groups
overall = overall.validateModel();
natural = natural.validateModel();


%% Set up boundary conditions and initial conditions
% We start by the boundary conditions. Here, we use standard routines from
% MRST, but must remember to also impose component specification at all
% boundaries
p = 50*barsa; s = []; z = [1, 0];     % p, T, z, s
bc = fluxside([], G, 'xmin', 1/day, 'sat', [0, 1]);    % Flux
bc = pside(bc, G, 'xmax', p, 'sat', [0, 1]);           % Standard bc
bc.components = repmat([0, 1], numel(bc.face), 1);     % Boundary z

% We first set the initial state from the prescribed parameters, validate
% the model (which will add data structures for any wells, etc), and then
% perform a flash to ensure that the corresponding initial guess exists at
% equilibrium conditions
state0 = eCPAinitCompositionalState(overall, p, T, s, z, false);
state0 = overall.validateState(state0);
initGuess = state0;
initGuess = overall.computeFlash(initGuess);

%% Adjust solver settings and construct nonlinear solver
% To improve the nonlinear solver we impose limits on changes during
% each nonlinear iteration. With natural variables, we limit the changes for
% phase mole fractions and saturation, whereas for the overall composition,
% we limit the overall mole fraction change.
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
[solNat1, reportsNat1] = nls.solveTimestep(state0, dt, natural, ...
                            'bc', bc, 'initialGuess', initGuess);
[solMole1, reportsMole1] = nls.solveTimestep(state0, dt, overall, ...
                            'bc', bc, 'initialGuess', initGuess);

%% Extract data from the output states
getStates = @(reports) cellfun(@(x) x.state, reports.StepReports{1}.NonlinearReport, 'UniformOutput', false);
% Get the outputs
getZ = @(states) [initGuess.components(1); cellfun(@(x) x.components(:, 1), states)];
getS = @(states) [initGuess.s(:, 1); cellfun(@(x) x.s(:, 1), states)];
getP = @(states) [initGuess.pressure(1); cellfun(@(x) x.pressure, states)] - p;
natStates = getStates(reportsNat1);
moleStates = getStates(reportsMole1);

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
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
