function setup = drogon_bo(varargin)
%Setup function for the Drogon model

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
        'The synthethic Drogon model from Equinor. See ', ...
        'https://webviz-subsurface-example.azurewebsites.net/', ...
        'drogon-conceptual-description ', ...
        'for details.'
    ];
    options = struct();
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
    % Check that OPM data is downloaded and registered as a module
    opm = mrstPath('opm-tests');
    assert(~isempty(opm), ['You must register ', ...
        'https://github.com/opm/opm-tests as a module using mrstPath in', ...
        'order to use this setup function']);
    % Make dataset path
    fn = fullfile(mrstPath('opm-tests'), 'drogon', 'model', 'DROGON_HIST.DATA');
    % Initialize initial state, model, and schedule
    [state0, model, schedule] = initEclipseProblemAD(fn, options.extra{:});
    % Normalize saturations
    state0.s = state0.s./sum(state0.s, 2);
    % Set defaulted threshold pressures
    model = model.validateModel();
    ppd = model.FlowDiscretization.PhasePotentialDifference;
    ppd = ppd.setThresholdPressuresFromState(model, state0);
    model.FlowDiscretization.PhasePotentialDifference = ppd;
    % Use simple wells
    for i = 1:numel(schedule.control)
        W = schedule.control(i).W;
        [W.isMS] = deal(false);
        W = rmfield(W, 'topo');
        schedule.control(i).W = W;
    end
    % Plotting
    plotOptions = { ...
        'View'              , [94, 27]                   , ...
        'Size'              , [750, 500]                 , ...
        'Box'               , true                       , ...
        'PlotBoxAspectRatio', [23.3129, 25.0413, 5.0000]   ...
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