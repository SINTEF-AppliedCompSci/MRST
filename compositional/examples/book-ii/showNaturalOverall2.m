%% Comparison of the Natural Variable/Overall Composition, Part 2
% This script contains a continuation of the showNaturalOverall script that
% compares the natural-variable and the overall-composition simulation. The
% setup of the reservoir is the same (a single cell), but this script is
% configured so that you can run different fluid mixtures instad of the
% simple brine-CO2 model considered in the first part.
%
% The setup can be configured through three variables:
%   - mixture_name gives the name of the model used for the fluid mixture.
%     See getBenchmarkMixture for possible values. Default: 'lumped_1'
%   - T_offset gives how much you want to offset the temperature from the
%     original value of 303.15 K. Default: 0
%   - p_offset gives hw much you want to offset the pressure from the
%     original value of 50 bar. Default: 0
mrstModule add compositional ad-core ad-props
mrstVerbose on

%% Construct the two simulation models
% Grid and petrophysical data in a single cell
G    = computeGeometry(cartGrid(1));      % Single cell 1 m^3 grid
rock = makeRock(G, 0.5*darcy, 0.5);

% Start with fluid properties from a standard black-oil, water-gas model
% and add extra mixture properties to turn this into a compositional model
f = initSimpleADIFluid('phases', 'wg', 'blackoil', false, 'rho', [1000, 700]);

if ~exist('mixture_name', 'var'), mixture_name = 'lumped_1'; end
if ~exist('T_offset', 'var'),     T_offset = 0;              end
if ~exist('p_offset', 'var'),     p_offset = 0;              end
[mixture, info] = getBenchmarkMixture(mixture_name);

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

%%
if ~exist('maxChange', 'var')
    maxChange = 0.1;
end
if isfinite(maxChange)
    % Limit changes - natural variables limits phase mole fractions and
    % saturation while overall composition limits the overall mole fraction
    % change
    overall.dzMaxAbs = maxChange;
    natural.dzMaxAbs = maxChange;
    natural.dsMaxAbs = maxChange;
end

%% Set up BC + initial conditions/initial guess
% We start by the boundary conditions. Here, we use standard routines from
% MRST, but must remember to also impose component specification at all
% boundaries
eos = overall.EOSModel;
p   = info.pressure + p_offset;
T   = info.T + T_offset;
bc  = fluxside([], G, 'xmin', 1/day, 'sat', [0, 1]);
bc  = pside(bc, G, 'xmax', p, 'sat', [0, 1]);
bc.components = repmat(info.injection, numel(bc.face), 1);

% We first set the initial state from the prescribed parameters, validate
% the model (which will add data structures for any wells, etc), and then
% perform a flash to ensure that the corresponding initial guess exists at
% equilibrium conditions
state0 = initCompositionalState(G, p, T, [1, 0], info.initial, eos);
state0 = overall.validateState(state0);

% Set up bad initial guess
initGuess = state0;
initGuess.components = zeros(1, eos.getNumberOfComponents());
initGuess.components(1) = 1;
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

% Set up nonlinear solver 
nls = NonLinearSolver('reportLevel', 3,'MaxIterations', 1000);

%% Compute a single time step
dt = 1*day;
[sn, nreports]  = nls.solveTimestep(state0, dt, natural, ...
                        'bc', bc, 'initialGuess', initGuess);
[sm, mreports] = nls.solveTimestep(state0, dt, overall, ...
                        'bc', bc, 'initialGuess', initGuess);
                    
%% Display runtime and number of iterations
figure
subplot(1,2,1), 
bar([nreports.Iterations; mreports.Iterations]); 
set(gca,'XTickLabel', {'Natural','Overall'}); title('Iterations');
subplot(1,2,2), 
bar([nreports.WallTime; mreports.WallTime]); 
set(gca,'XTickLabel', {'Natural','Overall'}); title('Wall time');

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
