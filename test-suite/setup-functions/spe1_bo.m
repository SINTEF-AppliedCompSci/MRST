function setup = spe1_bo(varargin)
% Setup function for the first SPE Comparative Solution Project (SPE 1)
%
% SYNOPSIS:
%   setup = spe1_bo('pn1', pv1, ...)
%   setup = spe1_bo(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup for the first SPE Comparative Solution Project (A.S. Odeh,
%   Comparison of solutions to a three-dimensional black-oil reservoir
%   simulation problem, J. Pet. Technol., 33(1):13-25, 1981, doi:
%   10.2118/9723-PA).
%
%   The reservoir model is represented on a 10x10x3 Cartesian grid and
%   describes gas injection into an undersaturated oil (i.e., an oil that
%   can dissolve more oil). Although this is essentially a two-phase
%   gas-oil case, it is represented as a three-phase problem with a water
%   phase having relative permeability virtually equal zero over the whole
%   saturation range. The reservoir is produced by a gas injector and a
%   bhp-controlled producer, placed on opposite ends of the diagonal.
%   Simulation of this case is discussed in more detail in section 11.8.1
%   of the first MRST textbook (Lie, Cambridge University Press, 2019, doi:
%   10.1017/9781108591416).
%
%   Configurable parameters: none.
%
% RETURNS:
%   setup - test case with the following fields: name, description,
%      options, state0, model, schedule, and plotOptions.
%      If the optional input fullSetup (see synopsis) is false, the
%      returned setup only contains name, description, and options.
%
% SEE ALSO:
%   TestCase, testcase_template, testSuiteTutorial, spe3_bo.


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
        ['SPE1 benchmark, see Odeh, A.S., "Comparison of Solutions to '            , ...
         'a Three-Dimensional Black-Oil Reservoir Simulation Problem.", '          , ...
         'J. Pet. Technol. 33 (1): 13-25. SPE-9723-PA. 1981, doi: 10.2118/9723-PA '];

    % Optional input arguments
    options = struct();
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end

    % Define module dependencies
    require ad-core ad-props ad-blackoil deckformat

    % Get initial state, model, and schedule
    pth = getDatasetPath('spe1');
    fn  = fullfile(pth, 'BENCH_SPE1.DATA');
    [state0, model, schedule] = initEclipseProblemAD(fn);
    
    % Set plot options
    plotOptions = {'PlotBoxAspectRatio', [1,1,0.25], ...
                   'View'              , [-15, 20] };
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end