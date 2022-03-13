function setup = qfs_compositional(varargin)
% Setup function for a quarter five-spot compositional model
%
% SYNOPSIS:
%   setup = qfs_compositional('pn1', pv1, ...)
%   setup = qfs_compositional(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup of a quarter of a five-spot case in which a supercritical carbon
%   dioxide containing some n-decane is injected into a homogeneous
%   reservoir initially filled by a mixture of n-decane, methane, and
%   carbon dioxide.
%
%   Confingurable parameters:
%      'ncells' - number of cells in x/y direction (default: 60)
%      'nsteps' - number of time steps             (default: 100)
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
    description = ...
        ['Quarter five-spot with CO2 injection into methane, n-decane, and CO2. '       , ...
         'Slightly modified from Klemetsdal et al, SPE RSC 2019, doi: 10.2118/193934-ms'];

    % Optional input arguments
    options = struct('ncells', 60, ... % Number of cells in x- and y-directions
                     'nsteps', 100);   % Number of timesteps
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end
    
    % Define module dependencies
    require ad-core ad-props compositional
    
    % Model
    G     = computeGeometry(cartGrid([options.ncells, options.ncells, 1], [1000, 1000, 10]));
    rock  = makeRock(G, 0.1*darcy, 0.25);
    fluid = initSimpleADIFluid('n'     , [2, 3]    ,...
                               'phases', 'OG'      , ...
                               'rho'   , [100, 100]);
    [cf, info] = getBenchmarkMixture('simple');
    model      = GenericOverallCompositionModel(G, rock, fluid, cf, 'water', false);
    
    % Schedule
    time  = 5*year;
    irate = 100*sum(model.operators.pv)/time;
    W = [];
    W = addWell(W, G, rock, 1, 'type', 'rate', 'val', irate, 'comp_i', [1, 0]);
    W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 25*barsa, 'comp_i', [1, 0]);
    for i = 1:numel(W)
        W(i).components = info.injection;
    end
    dt       = rampupTimesteps(time, time/options.nsteps);
    schedule = simpleSchedule(dt, 'W', W);
    
    % Initial state
    state0 = initCompositionalState(G, 50*barsa, info.temp, [1, 0], info.initial, model.EOSModel);
    
    % Plotting
    plotOptions = {'PlotBoxAspectRatio', [1,1,1]       , ...
                   'View'              , [0, 90]       , ...
                   'Projection'        , 'orthographic'};
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end