function setup = testcase_template(varargin)
%Template for setup functions from the test-suite module.
%
% SYNOPSIS:
%   setup = testcase_template('pn1', pv1, ...)
%   setup = testcase_template(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Template for all setup functions in the test suite. Function names
%   should be lower-case, with underscore to distinguish words. If
%   possible, the name should be on the format 'simple_description_phases'
%   or 'simple_description_compositional', e.g., 'my_testcase_wog'. May be
%   used as a standalone test-case definition, or to construct an instance
%   of the TestCase class as testcase = TestCase('my_testcase_wog');
%
% OPTIONAL PARAMETERS:
%   The setup function may take in any number of optional input arguments
%   to facilitate different configurations of the same test case, e.g.,
%   number of cells, number of timesteps, and so on:
%   Example: 
%     my_testcase_wog('ncells', 2020, 'nsteps', 5);
%   Suggested standard names:
%     'ncells' - Number of cells, e.g. in each axial direction, or total
%     'nsteps' - Number of timesteps
%     'time'   - Total simulation time
%     'dt'     - Target timestep length
%
%   If called with a boolean first input argument, this is interpreted as a
%   flag indicating whether the function shall return options and
%   description only. This is usesful for getting a test-case description
%   without having to construct the entire simulation model and is, for
%   example, used by the TestCase class to check whether the test case is
%   already stored to disk.
%
% RETURNS:
%   setup - Test case struct with the following fields:
%           - name - Name of the test case. By default, the name is the
%                 same as the function itself, but can also be chosen to
%                 explain parameter settings.
%
%           - description - One-line test case description
%
%           - options - A struct of the optional input parameters, with
%                 defaults for all arguments that were not passed as
%                 optional parameters. Returned for convenient access to
%                 the test case configuration.
%
%           - state0, model, schedule - Initial state, model, and
%                 simulation schedule that can be passed to
%                 simulateScheduleAD or packSimulationProblem
%
%           - plotOptions - Cell array on the form {'pn1', pv1, ...}
%                 with arguments that can be used in any of the following
%                 ways
%                   - set(myAxis, 'pn1, vn1, ...)
%                   - figure('pn1', vn1, ...)
%                   - plotToolbar(G, state, 'pn1', vn1, ...)
%                 In addition to the standard optional parameters of
%                 `figure`, {'Size', [width, height]} can also be provided,
%                 which the TestCase class interprets as
%                 [pos(1:2), [width, height]], where
%                 pos = get(0, 'DefaultFigurePosition')
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

    % A test case setup function consists of three parts: definition of
    % description and options, definition of module dependencies, setting
    % up the initial state, model schedule and options for plotting, and
    % packing the test case setup.
    
    % Step 1: Test case description and options
    %---------------------------------------------------------------------%
    % Always include this block (with obvious modifications)
    %---------------------------------------------------------------------%
    % Each test case must start with the description and options, followed
    % by an nargout check that returns if we only asked for the description
    % and options.
    description = 'Test case template. Not meant for direct use';
    % Each test case can have any number of optional input arguments, the
    % setup returned must include a (possibly empty) options struct
    options = struct('ncells', 2020, ...
                     'nstep' , 5   );
    % Process optinal input arguments
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
                                        options, description, varargin{:});
    if ~fullSetup, return; end
    %---------------------------------------------------------------------%
    
    % Step 2: Define any module dependencies for the test case and set up
    %---------------------------------------------------------------------%
    % initial state, model and schedule.
    % The following dependencies are typically always needed
    require ad-core ad-props ad-blackoil
    % Define initial state, model and schedule
    [state0, model, schedule] = deal([]);
    % plotOptions are only used by TestCase. In case of empty
    % plotOptions, TestCase will attempt to set reasonable defaults
    plotOptions = {};
    %---------------------------------------------------------------------%
    
    % Step 3: Pack test case setup
    %---------------------------------------------------------------------%
    % Always include this block (with obvious modifications)
    %---------------------------------------------------------------------%
    % Pack setup. This is done using the packTestCaseSetup, which takes in
    % a test case name, in addition to the optional input arguments listed.
    % All optional arguments not passed to the function will default to
    % empty fieilds in the setup.
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
    %---------------------------------------------------------------------%
    
end