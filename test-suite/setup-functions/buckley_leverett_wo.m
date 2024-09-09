function setup = buckley_leverett_wo(varargin)
% Setup function for standard 1D Buckley-Leverett problem
%
% SYNOPSIS:
%   setup = buckley_leverett_wo('pn1', pv1, ...)
%   setup = buckley_leverett_wo(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup of a standard 1D Buckley-Leveret case an immisicble two-phase
%   fluid. The configurable parameters include:
%      'ncells'    - number of cells in x/y direction (default: 100)
%      'nkr'       - Brooks-Corey relperm exponent    (default: 2)
%      'cfl'       - CFL number relating dt and dx    (default: 1)
%      'realistic' - Use typical reservoir values for parameters. By
%      default, this routine uses units that are 1 or close to 1. This is
%      very different than typical reservoir values in SI. Enabling
%      realistic adds weak compressibility and makes the case closer to one
%      that can be run in a standard reservoir simulator, at the cost of
%      making the parameters less legible from a mathematical point of
%      view.
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
    % One-line description
    description = 'Two-phase Buckley Leverett displacement';
    % Optional input arguments
    options = struct('ncells',    100, ... % Number of cells
                     'nkr'   ,    2  , ... % Brooks-Corey relperm exponent
                     'realistic', false, ...
                     'cfl'   ,    1  );    % CFL number
    [options, fullSetup, setup] = processTestCaseInput(mfilename, options, description, varargin{:});
    if ~fullSetup, return; end
    % Define module dependencies
    require ad-core ad-props ad-blackoil
    % Model
    if options.realistic
        pressure = 500*barsa;
        perm = 0.1*darcy;
        poro = 0.1;
        L = 1000*meter;
        c = [1e-6, 1e-5]/barsa;
        WI = -1;
        mu = 0.1*centi*poise;
        tscale = day;
    else
        [perm, poro, pressure, mu, tscale] = deal(1);
        L = options.ncells;
        c = [0, 0];
        WI = 9999;
    end
    G     = computeGeometry(cartGrid([options.ncells,1], [L, 1]));
    rock  = makeRock(G, perm, poro);
    fluid = initSimpleADIFluid('phases', 'WO', 'n', [1,1]*options.nkr, 'mu', [1,1]*mu, 'rho', [1,1], 'c', c, 'pRef', pressure);
    model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
    % Set up wells
    W = [];
    time = tscale*options.ncells/options.cfl;
    rate = sum(poreVolume(G, rock))/time;
    W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate, ...
                                         'WI', WI, 'comp_i', [1,0], 'sign', 1);
    W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', pressure/2, ...
                                         'WI', WI, 'comp_i', [1,0], 'sign', -1);
    % Make schedule
    schedule = simpleSchedule(rampupTimesteps(time, 1*tscale, 0), 'W', W);
    % Initial state
    state0 = initResSol(G, pressure, [0,1]);
    % Plotting
    plotOptions = {'plot1d'            , true              , ...
                   'lockCaxis'         , true              , ...
                   'Size'              , [800, 350]        , ...
                   'PlotBoxAspectRatio', [2.5,1,1]         , ...
                   'XLim'              , [0,options.ncells], ...
                   'YLim'              , [0,1]             };
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end