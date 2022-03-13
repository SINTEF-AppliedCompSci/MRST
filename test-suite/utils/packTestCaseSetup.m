function setup = packTestCaseSetup(name, varargin)
%Utility function for packing test case setup.
%
% SYNOPSIS:
%   setup = packTestCaseSetup(name)
%   setup = packTestCaseSetup(name, 'pn1', pv1, ...)
%
% REQUIRED PARAMETERS:
%   name - name of test case as string
%
% OPTIONAL PARAMETERS:
%   description             - One-line test case description
%   options                 - Test case options struct
%   state0, model, schedule - Test case setup
%   plotOptions             - Options for plotting
%   extra                   - Field for anything that does not fit into the
%                             categories above (e.g., reference results)
%
% RETURNS:
%   setup - Test case struct with all parameters described above
%
% SEE ALSO:
%   `testcase_template`, `TestCase`, `processTestCaseInput`

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

    setup = struct('name'       , name, ...
                   'description', []  , ...
                   'options'    , []  , ...
                   'state0'     , []  , ...
                   'model'      , []  , ...
                   'schedule'   , []  , ...
                   'plotOptions', {{}}, ...
                   'extra'      , []  );
    setup = merge_options(setup, varargin{:});
    
end