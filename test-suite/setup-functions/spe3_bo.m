function setup = spe3_bo(varargin)
% Setup function for the first SPE Comparative Solution Project (SPE 1)
%
% SYNOPSIS:
%   setup = spe3_bo('pn1', pv1, ...)
%   setup = spe3_bo(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup for the third SPE Comparative Solution Project (D.E. Kenyon and
%   G.A. Behie, Third SPE Comparative Solution Project: Gas cycling of
%   retrograde condensate reservoirs, JPT, August 1987 (981)', doi:
%   10.2118/12278-PA).
%
%   The model is a live-gas black-oil version of the SPE 3 benchmark, which
%   originally described a compositional model for a condensate reservoir.
%   The reservoir model has 9x9x4 cells, describing four layers of
%   different permeability and thickness.
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
%   TestCase, testcase_template, testSuiteTutorial, spe1_bo


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
        ['SPE3 benchmark, see Kenyon, D.E. and Behie, G.A., '       , ...
         '"Third SPE Comparative Solution Project: Gas Cycling of ' , ...
         'Retrograde Condensate Reservoirs," JPT, August 1987 (981)', ...
         'doi:10.2118/12278-PA'                                     ];

    % Optional input arguments
    options = struct();
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end
    
    % Define module dependencies
    require ad-core ad-props ad-blackoil deckformat
    
    % Get initial state, model, and schedule
    pth = getDatasetPath('spe3');
    fn  = fullfile(pth, 'BENCH_SPE3.DATA');
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