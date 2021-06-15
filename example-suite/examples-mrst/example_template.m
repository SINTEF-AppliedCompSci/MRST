function [description, options, state0, model, schedule, plotOptions] = example_template(varargin)
%Example template for the example suite. 
%
% SYNOPSIS:
%   [description, options, state0, model, schedule, plotOptions] = example_template('pn1', pv1, ...)
%
% DESCRIPTION:
%   Template for all examples in the example suite. Function names should
%   be lower-case, with underscore to distinguish words. If possible, the
%   name should be on the format 'simple_description_phases' or
%   'simple_description_compositional', e.g., 'my_example_wog'. May be used
%   as a stand-alone example definition, or to construct an instance of
%   `MRSTExample` as example = MRSTExample('my_example_wog');
%
% OPTIONAL PARAMETERS:
%   The example function may take in any number of optional input
%   arguments to facilitate different configurations of the same example,
%   e.g., number of cells, number of timesteps, and so on:
%   Example: my_example_wog('ncells', 2020, 'nsteps', 5);
%   Suggested standard names:
%   
%   'ncells' - Number of cells, e.g. in each axial direction, or total
%   'nsteps' - Number of timesteps
%   'time'   - Total simulation time
%   'dt'     - Target timestep length
%
% RETURNS:
%   description - One-line example description, displayed in list-examples,
%                 and the only input argument if the function is called as
%                 description = my_example_wog()
%
%   options     - A struct of the optional input parameters, with defaults
%                 for all arguments that were not passed as optional
%                 parameters. Returned for convenient access to the example
%                 configuration.
%
%   state0, model, schedule - Initial state, model, and simulation schedule
%                             that can be passed to `simulateScheduleAD`
%
%   plotOptions - Cell array on the form {'pn1', pv1, ...} with arguments
%                 that can be used in any of the following ways
%                   - set(myAxis, 'pn1, vn1, ...)
%                   - figure('pn1', vn1, ...)
%                   - plotToolbar(G, state, 'pn1', vn1, ...)
%                 In addition to the standard optional parameters of
%                 `figure`, {'Size', [width, height]} can also be provided,
%                 which `MRSTExample` interprets as
%                 [pos(1:2), [width, height]], where
%                 pos = get(0, 'DefaultFigurePosition')
%
% SEE ALSO:
%   `MRSTExample`, `listExamples`, `exampleSuiteTutorial`

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
    % Each example must start with the description and options, followed by
    % an nargout check that returns if we only asked for the description
    % and options
    description = 'Example template. Not meant for direct use';
    % Each example can have any number of optional input arguments, and
    % must return a (possibly empy) options struct
    options = struct('ncells', 2020, 'nstep', 5);
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
    % Define any module dependencies for the example. The following are
    % typically always needed
    require ad-core ad-props ad-blackoil
    % Define initial state, model and schedule
    [state0, model, schedule] = deal([]);
    % plotOptions are only by MRSTExample. In case of empty plotOptions,
    % MRSTExample will attempt to set reasonable defaults
    plotOptions = {};
end