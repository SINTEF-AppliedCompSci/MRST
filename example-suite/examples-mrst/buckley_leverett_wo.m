function [description, options, state0, model, schedule, plotOptions] = buckley_leverett_wo(varargin)
%Example from the example suite, see description below.
%
% SEE ALSO:
%   `MRSTExample`, `example_template`, `exampleSuiteTutorial`.

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
    % One-line description
    description = 'Two-phase Buckley Leverett displacement';
    % Optional input arguments
    options = struct('ncells', 100, ... % Number of cells
                     'nkr'   , 2  , ... % Brooks-Corey relperm exponent
                     'cfl'   , 1  );    % CFL number
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
    % Define module dependencies
    require ad-core ad-props ad-blackoil
    % Model
    G     = computeGeometry(cartGrid([options.ncells,1]));
    rock  = makeRock(G, 1, 1);
    fluid = initSimpleADIFluid('phases', 'WO', 'n', [1,1]*options.nkr, 'mu', [1,1], 'rho', [1,1]);
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    % Set up wells
    W = [];
    time = options.ncells/options.cfl;
    rate = sum(poreVolume(G, rock))/time;
    W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate, 'WI', 9999, 'comp_i', [1,0]);
    W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 0.5, 'WI', 9999, 'comp_i', [1,0]);
    % Make schedule
    schedule = simpleSchedule(rampupTimesteps(time, 1, 0), 'W', W);
    % Initial state
    state0 = initResSol(G, 1, [0,1]);
    % Plotting
    plotOptions = {'plot1d'            , true              , ...
                   'lockCaxis'         , true              , ...
                   'Size'              , [800, 350]        , ...
                   'PlotBoxAspectRatio', [2.5,1,1]         , ...
                   'XLim'              , [0,options.ncells], ...
                   'YLim'              , [0,1]             };
end