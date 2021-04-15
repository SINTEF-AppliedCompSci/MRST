function [description, options, state0, model, schedule, plotOptions] = qfs_wo(varargin)
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
    description = 'Quarter five-spot example with two-phase fluid on Cartesian grid';
    % Optional input arguments
    options = struct('ncells', 50    , ... % Number of cells in x- and y-directions
                     'time'  , 2*year, ... % Total injection time
                     'dt'    , 30*day, ... % Timestep length
                     'pvi'   , 1     , ... % Pore volumes injected
                     'nkr'   , 2     );    % Brooks-Corey relperm exponent
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
    % Define module dependencies
    require ad-core ad-props ad-blackoil
    % Model
    G     = computeGeometry(cartGrid([1,1]*options.ncells, [1000, 1000]*meter));
    rock  = makeRock(G, 100*milli*darcy, 0.4);
    fluid = initSimpleADIFluid('phases', 'WO', 'n', [1,1]*options.nkr, 'mu', [1,1]*centi*poise, 'rho', [1,1]);
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    % Wells
    rate = options.pvi*sum(poreVolume(G, rock))/options.time;
    bhp  = 50*barsa;
    W = [];
    W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate, 'comp_i', [1,0]);
    W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', bhp, 'comp_i', [1,0]);
    % Schedule
    schedule = simpleSchedule(rampupTimesteps(options.time, options.dt), 'W', W);
    % Initial state
    state0      = initResSol(G, 10*barsa, [0,1]);
    % Default plotting
    plotOptions = {};
end