function setup = gravity_segregation(varargin)
% Setup function for gravity segregation test cases
%
% SYNOPSIS:
%   setup = gravity_segregation('pn1', pv1, ...)
%   setup = gravity_segregation(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup of for two types of gravity segregation test cases. The first
%   type consists of 1x1xN cells with horzontal fluid contacts separating
%   heavier and lighter fluid. For a three-phase system (WOG), the heaviest
%   fluid fills the upper 1/3 of the column, the lightest fluid the bottom
%   1/3, and the intermediate phase the remaining 1/3. For a two-phase
%   system, the heavier fluid fills the upper 1/2 of the volume.
%
%   The second type consists of a Nx1xN reservoir with vertical fluid
%   contacts separating heavier and lighter phases, ordered similarly as
%   for the first case.
%
%   The configurable parameters include:
%      'phases'  - string specifying phase present   (default: 'WOG')
%      'ncells'  - number of cells in each direction (default: 21)
%      'nsteps'  - number of time steps              (default: 100)
%      'type'    - horizontal/vertical fluid contact (default: 'horizontal')
%      'useRampup' - ramp up time steps              (default:  true)
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
    description ...
        = ['Gravity segregation example with initially horizontal or ', ...
           'vertical fluid contact and two- or three-phase immiscible fluids'];
    
    % Optional input arguments
    options = struct('phases'   , 'WOG'       , ...
                     'ncells'   , 21          , ...
                     'nsteps'   , 100         , ...
                     'type'     , 'horizontal', ...
                     'useRampup', true        );
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    
    % Set parameters depending on fluid contact-type
    nph = numel(options.phases);
    options.ncells = options.ncells - rem(options.ncells, nph);
    if ~fullSetup, setup.options = options; return; end
    
    % Define module dependencies
    require ad-core ad-props ad-blackoil
    L = 100*meter;
    switch options.type
        case 'vertical'
            cartDims   = [options.ncells, 1, options.ncells];
            physDims   = [L, 1, L];
            pbar       = [1,1,1];
            sz         = [500, 500];
            view       = [0,0];
        case 'horizontal'
            cartDims   = [1, 1, options.ncells];
            physDims   = [1, 1, L];
            pbar       = [0.2,0.2,1];
            sz         = [300, 500];
            view       = [120,25];
    end
    
    % Make model
    gravity reset on
    G    = cartGrid(cartDims, physDims); % Cartesian grid
    G    = computeGeometry(G);
    rock = makeRock(G, 1*darcy, 0.3);    % Homogeneous rock

    % Simple, incompresible fluid
    phases    = 'WOG';
    fluidArgs = {'phases', phases                            , ...
                 'mu'    , [1,5,1]*centi*poise               , ...
                 'rho'   , [1500, 1000, 500]*kilogram/meter^3, ...
                 'n'     , [2,2,2]                           };
    active = ismember(phases, options.phases);
    fluidArgs(2:2:end) = cellfun(@(arg) ...
        arg(active), fluidArgs(2:2:end), 'UniformOutput', false);
    fluid = initSimpleADIFluid(fluidArgs{:}, 'cr', 1e-8);
    
    % Black-oil model
    model = GenericBlackOilModel(G, rock, fluid, 'water', active(1), ...
                                                 'oil'  , active(2), ...
                                                 'gas'  , active(3));

    % Make schedule of five years with options.nsteps timesteps
    time = 5*year;
    dt   = rampupTimesteps(time, time/options.nsteps, 5*options.useRampup);
    schedule = simpleSchedule(dt);
    [ii, ~, kk] = gridLogicalIndices(G);
    switch options.type
        case 'vertical'
            ll = ii;
        case 'horizontal'
            ll = kk;
    end
    
    % Initial state
    state0 = initResSol(G, 1*atm, zeros(1,nph));
    d = 1/nph;
    for i = 1:nph
        cells = ll > max(ll)*d*(i-1) & ll <= max(ll)*d*i;
        sat = zeros(1,nph);
        sat(i) = 1;
        state0.s(cells,:) = repmat(sat, nnz(cells), 1);
    end
    
    % Plotting
    plotOptions = {'PlotBoxAspectRatio', pbar         , ...
                   'Projection'        , 'perspective', ...
                   'View'              , view         , ...
                   'Size'              , sz           };
    
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end