function setup = qfs_wo(varargin)
% Setup function for a quarter five-spot water-oil model
%
% SYNOPSIS:
%   setup = qfs_wo('pn1', pv1, ...)
%   setup = qfs_wo(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup of a standard quarter of a five-spot case for an immisicble
%   two-phase fluid. The configurable parameters include:
%      'ncells' - number of cells in x/y direction (default: 50)
%      'time'   - total simulation horizon         (default: 2 years)
%      'dt'     - uniform time step size           (default: 30 days)
%      'pvi'    - total poree volumes injected     (default: 1)
%      'nkr'    - Brooks-Corey relperm exponent    (default: 2)
%
% RETURNS:
%   setup - test case with the following fields: name, description,
%      options, state0, model, schedule, and plotOptions.
%      If the optional input fullSetup (see synopsis) is false, the
%      returned setup only contains name, description, and options.
%
% SEE ALSO:
%   TestCase, testcase_template, testSuiteTutorial.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end
    
    % Define module dependencies
    require ad-core ad-props ad-blackoil

    % Cartesian grid covering a 1000 x 1000 m domain
    G     = computeGeometry(cartGrid([1,1]*options.ncells, [1000, 1000]*meter));
    rock  = makeRock(G, 100*milli*darcy, 0.4); % Rock with 100 md perm and 0.4 poro
    fluid = initSimpleADIFluid('phases', 'WO'             , ... % Water and oil
                               'n'     , [1,1]*options.nkr, ... % Relperm exponents
                               'mu'    , [1,1]*centi*poise, ... % Viscosity
                               'rho'   , [1,1]            );    % Density (no gravity)

    % Construct two-phase model with oil and water
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
    state0 = initResSol(G, 10*barsa, [0,1]);
    
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   );
end