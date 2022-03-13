function [options, fullSetup, setup] = processTestCaseInput(name, options, description, varargin)
%Utility function for processing test case input.
%
% SYNOPSIS:
%   setup = processTestCaseInput(name, options, description)
%   setup = processTestCaseInput(fullSetup, name, options, description)
%   setup = processTestCaseInput(... , 'pn1', pv1, ...)
%
% REQUIRED PARAMETERS:
%   name       - Name of test case as string
%   options    - Test case default options struct
%   decription - One-line test case description
%
% OPTIONAL PARAMETERS:
%   fullSetup - Flag indicating if we intend to set up the entire test
%               case or not. If false, the third return variable is a test
%               case setup with only the name, description and options.
%               Default value is true.
%   Test case options - The function processes optional parameters for the
%               test case. Parameters that are not part of the options
%               struct are stored as a cell array on the form {'pn1', vn1,
%               ...} to options.extra to facilitate passing these to
%               functions called from inside the test case setup function.
%
% RETURNS:
%   options   - Test case options struct with user-defined options
%   fullSetup - Flag indicating if we should build the full test case
%   setup     - If fullSetup is false, this is a simplified test case setup
%               struct with only name, description and user-defined
%               options. If fullSetup is true, this is empty.
%   
%
% SEE ALSO:
%   `testcase_template`, `TestCase`, `packTestCaseSetup`

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

    % Determine if we should build the full setup
    fullSetup = true;
    if nargin > 3 && islogical(varargin{1})
        fullSetup = varargin{1}; varargin = varargin(2:end);
    end
    % Process options
    [options, extra] = merge_options(options, varargin{:});
    % Set any extra input to options struct
    options.extra = extra;
    if fullSetup
        % We will do the full setup later - return an empty setup
        setup = [];
    else
        % Pack simplified setup
        setup = packTestCaseSetup(name,                       ...
                                  'description', description, ...
                                  'options'    , options    );
    end
    % Issue warning to tell the user about optional arguments passed that
    % are not part of default test case options
    if ~isempty(extra)
        pl = ''; if numel(extra(1:2:end)) > 1, pl = 's'; end
        warning(['Option%s %s \n not part of default test case ', ...
                 'options. Storing to options.extra'],            ...
                 pl, sprintf('\n * %s', extra{1 : 2 : end}))
    end
    
end