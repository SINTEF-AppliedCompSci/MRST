function params = addParameter(params, SimulatorSetup, varargin)
% Simple utility function for adding a new parameter to a list of
% parameters.
% See  ModelParameter.m for more information about the parameters class.
%

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

p = ModelParameter(SimulatorSetup, varargin{:});
if isempty(params)
    params = {p};
elseif iscell(params)
    params = [params, {p}];
elseif isa(params, class(p))
    params = [{params}, {p}];
else
    error('Unknown format of input ''params''');
end
end
