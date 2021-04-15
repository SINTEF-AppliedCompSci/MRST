function opt = merge_options_relaxed(opt, varargin)
% A less general version of merge_options focused on specific choices:
%   - Arguments must match the names of fields exactly
%   - No type checking
%   - Errors for unsupported fields
%
% INTENTIONALLY UNDERDOCUMENTED, SUBJECT TO CHANGE.

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

    if nargin == 2 && iscell(varargin{1})
        % Can pass in varargin cell array directly
        varargin = varargin{1};
        n = numel(varargin);
    else
        % varargin is a normal list
        n = nargin - 1;
    end
    assert(mod(n, 2) == 0);
    for i = 1:2:(n-1)
        f = varargin{i};
        v = varargin{i+1};
        assert(~isstruct(f) || isfield(opt, f));
        opt.(f) = v;
    end
end
