%% Introduction to Vapor-Liquid Equilibrium Calculations
% This script contains examples from Section 8.3 of the second MRST book:
% Advanced Modelling with the MATLAB Reservoir Simulation Toolbox (MRST),
% Cambridge University Press, 2021.
mrstModule add compositional ad-core

%% §8.3.1: Show Rachford-Rice
% In the first example, we demonstrate how to use the Rachford-Rice
% algorithm to determine the vapor-liquid equilibrium for a two-component
% fluid system consisting of heavy and light molecules. The system has
% equilibrium constants Kh= 1/10,Kl= 10. Under two-phase conditions, the
% heavy component will thus mostly be present in the liquid phase, with one
% mole in the vapor phase for every ten moles in the liquid. The light
% component reverses the situation, with one out of eleven molecules ending
% up in the liquid phase.
%
% We span the range of possible mole fractions for the light component and
% solve the objective function to determine the corresponding vapor mole
% fraction at equilibrium.
n = 50;                         % Number of samples
K = [0.1, 10];                  % First component is heavy, second is light
z_light = linspace(0, 1, n)';   % Go from 0 -> 1
L = solveRachfordRiceVLE([], K, [1 - z_light, z_light]); % No initial guess
figure
plot(z_light, 1 - L, 'o-', 'linewidth', 1.5)
xlabel('z_l (light component)'); ylabel('V');

%%
% We plot the objective functions for K-value mass balance as function of
% vapor mole fraction at zl=zh= 0.5.
lw = 1.5;
z = [0.5, 0.5];
V = linspace(-1, 2, n)';
R = sum(((K - 1).*z)./(1 + V.*(K - 1)), 2);

R_1 = sum((K).*z./(1 + V.*(K - 1)), 2);
R_2 = sum(z./(1 + V.*(K - 1)), 2);

clf;
plot(V, R_1, '--', 'linewidth', lw); hold on; grid on
plot(V, R_2, '--', 'linewidth', lw); hold on; grid on
plot(V, R, '-k', 'linewidth', 1.8); hold on; grid on

plot(0.5, 0, 'ko', 'markerfacecolor', 'r')
xlabel('Vapor mole fraction')
ylabel('Objective')
legend('O_1', 'O_2', 'O_{RR}', 'Solution')


%% §8.3.2: Updating the Thermodynamic Equilibrium
% To use the Rachford-Rice method to determine the vapor-liquid
% equilibrium, we must know the K-values and suitable relationships for the
% densities. In the next example, we show how one can use an equation of
% state (Peng-Robinson) to determine the equilibrium state instead. We
% start by constructing the necessary input parameters to the PR EoS for a
% system consisting of three named components. These properties are
% generated from CoolProp and we display them using a custom-made disp
% function.
names = {'CarbonDioxide', 'Methane', 'n-Decane'};
mixture = TableCompositionalMixture(names);
disp(mixture)

peng_robinson = EquationOfStateModel([], mixture, 'Peng-Robinson');

%%
% Having constructed the equation of state, we use it to compute a phase
% diagram of the liquid fraction for the three-component system. The points
% in the diagram are computed by solving a standalone flash for the
% prescribed PR EoS
%e = 1e-3;
%single = [1, 1, 1];

% Generate compositions that span the phase diagram
x = linspace(0, 1, 100);
[X, Y] = meshgrid(x, x);
Z = 1 - X - Y;
c = [reshape(X, [], 1), reshape(Y, [], 1), reshape(Z, [], 1)];
c = bsxfun(@rdivide, c, sum(c, 2));

% Perform the standalone flash calculation to determine liquid fraction 
p = 25*barsa;
T = (273.15 + 30)*Kelvin;
[L, x, y] = standaloneFlash(p, T, c, peng_robinson);                       %#ok<ASGLU>

% Filter out any negative Z-values and plot ternary diagram
L(Z<-1e-8) = nan;
figure;
short_names  = {'CO_2', 'C_1H_4', 'C_{10}H_{22}'}; % Make plot nicer
[mapx, mapy] = ternaryAxis('names', short_names, 'isox', false, 'isoy', false);
xt = mapx(X, Y, Z);
yt = mapy(X, Y, Z);
L  = reshape(L, size(xt));
L(L==1 | L==0) = nan;
contourf(xt, yt, L, 10)
colorbar('horiz'); caxis([0, 1]); axis equal tight; ylim([0, 1])
cmap = interp1([1; 100], [0, 0, 1; 0.9, 0.9, 1], linspace(1, 100, 10));
colormap(cmap)

%%
% Compute similar phase diagrams for the fluid model from the SPE 5
% benchmark case.

% Load the fluid parameters and display them
[spe5, info] = getBenchmarkMixture('spe5');
disp(spe5)

% Construct the equation of state
eos = EquationOfStateModel([], spe5);
eos.extraOutput = 1;
eos.verbose = true;

% Perform the standalone flash
ns = 500;
z  = info.initial;
p  = linspace(1, 180, ns)*barsa;
t  = 273.15 + linspace(1, 400, ns);
[T, P] = meshgrid(t, p);
[L, x, y, ~, ~, ~, ~, reports] = standaloneFlash(P, T, z, eos);

%%
% Visualize the (T,p) phase diagram
clf;
L = reshape(L, ns, ns);
L(L==1 | L==0) = nan;
contourf(T - 273.15, P/barsa, L);
shading flat; view(0, 90);
xlabel('T [Celsius]')
ylabel('Pressure [Bar]');
caxis([0, 1]); colormap(cmap); colorbar

%%
% Visualize the number of iterations in the SSI method
its = 0;
for i = 1:numel(reports)
    its = its + double(reports{i}.ActiveFlag);
end
its(its == 0) = nan;
figure;
[C, h] = contourf(T - 273.15, P/barsa, reshape(its, ns, ns));
colorbar()
shading flat; view(0, 90);
xlabel('T [Celsius]');
ylabel('Pressure [Bar]');

%%
% Accelerate the SSI algorithm by switching to Newton after 5 iterations
% and plot a convergence of how fast the two methods converge
eos.maxSSI = 5;
eos.AutoDiffBackend = DiagonalAutoDiffBackend('rowMajor', true, 'usemex', true);
[L_n, ~, ~, ~, ~, ~, ~, reports_newton] = standaloneFlash(P, T, z, eos);


% Plot the convergence histories of the two flash calculations
ssi = cellfun(@(x) x.ActiveCells, reports);
newt = [cellfun(@(x) x.ActiveCells, reports_newton); 0];
figure; hold on
% yyaxis left
plot(ssi, '.')
plot(newt, '.', 'MarkerSize', 10)
set(gca, 'YScale', 'log')
xlabel('Iterations')
ylabel('Unconverged values')
legend('SSI', 'SSI+Newton');

%%
% Plot the time consumption for the flash and stability calculations
t_ssi       = sum(cellfun(@(x) x.TotalTime, reports));
t_newt      = sum(cellfun(@(x) x.TotalTime, reports_newton));
t_stability = sum(cellfun(@(x) x.StabilityTime, reports));

clf;
bar([t_stability, t_stability; t_ssi, t_newt]')
set(gca, 'XTickLabel', {'SSI', 'SSI+Newton'});
legend('Stability', 'Flash')
ylabel('Total time [s]')


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