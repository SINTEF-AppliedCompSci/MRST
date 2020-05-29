function xv = value(x)
%Remove AD state and compact 1 by n cell arrays to matrices
%
% SYNOPSIS:
%    v = value(V);
%
% DESCRIPTION:
%   Removes AD variables, and if the input is a cell array with a single
%   row and multiple columns, will combine them into a single matrix.
%
% REQUIRED PARAMETERS:
%   v   - Value to be converted.
%
% RETURNS:
%   v    - Value with no AD status.
%
% SEE ALSO:
%   double2ADI, double2GenericAD, ADI, GenericAD

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    if isnumeric(x) || islogical(x)
        xv = x;
    elseif iscell(x)
        sz = size(x);
        xv = cellfun(@value, x, 'UniformOutput', false);
        if sz(1) == 1
            % Cell arrays are converted to matrices
            xv = [xv{:}];
        end
    elseif isstruct(x)
        fn = fieldnames(x);
        for i = 1:numel(fn)
            f = fn{i};
            for j = 1:numel(x)
                x(j).(f) = value(x(j).(f));
            end
        end
        xv = x;
    else
        xv = x;
    end
end