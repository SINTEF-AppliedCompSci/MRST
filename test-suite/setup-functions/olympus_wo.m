function setup = olympus_wo(varargin)
%Setup function for the Olympus benchmark

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

    % Step 1: Test case description and options
    %---------------------------------------------------------------------%
    description = [ ...
        'The Olympus field optimization challenge, with 50 '        , ...
        'realizations of the same field model. '                    , ...
        'See Fonseca et al., Comput. Geosci. 24, 1933–1941 (2020), ', ...
        'doi: 10.1007/s10596-020-10003-4'
    ];
    options = struct( ...
        'realization', 1 ...
    );
    % Process optinal input arguments
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
                                        options, description, varargin{:});
    if ~fullSetup, return; end
    %---------------------------------------------------------------------%
    
    % Step 2: Define any module dependencies for the test case and set up
    %---------------------------------------------------------------------%
    % Module dependencies
    require ad-core ad-props ad-blackoil
    gravity reset on
    
    % Find dataset path
    assert(options.realization >= 1 && options.realization <= 50, ...
        'Realization must be between 1 and 50');
    pth  = fullfile(getDatasetPath('olympus'), 'ECLIPSE');
    name = sprintf('OLYMPUS_%d', options.realization);
    fn   = fullfile(pth, name, sprintf('%s.DATA', name));
    % Initialize initial state, model, and schedule
    [state0, model, schedule] = initEclipseProblemAD(fn, options.extra{:});
    
    % Plotting
    plotOptions = { ...
        'View'              , [15, 30]                   , ...
        'Size'              , [750, 500]                 , ...
        'Box'               , true                       , ...
        'PlotBoxAspectRatio', [39.8184, 19.8526, 10.0000]  ...
    };
    %---------------------------------------------------------------------%
    
    % Step 3: Pack test case setup
    %---------------------------------------------------------------------%
    setup = packTestCaseSetup(mfilename, ...
        'description', description, ...
        'options'    , options    , ...
        'state0'     , state0     , ...
        'model'      , model      , ...
        'schedule'   , schedule   , ...
        'plotOptions', plotOptions  ...
    );
    %---------------------------------------------------------------------%
    
end