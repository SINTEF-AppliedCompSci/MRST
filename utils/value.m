function xv = value(x, flatten)
%Remove AD state and optionally compact 1 by n cell arrays to matrices
%
% SYNOPSIS:
%    v = value(V);        % Flatten cell arrays with single row, multi col
%    v = value(V, false); % Don't flatten cell arrays
%
% DESCRIPTION:
%   Removes AD variables, and if the input is a cell array with a single
%   row and multiple columns, will combine them into a single matrix.
%
% REQUIRED PARAMETERS:
%   v       - Value to be converted.
%
% OPTIONAL PARAMETERS:
%   flatten - If true, any cell array inputs with a single row will be
%             converted to a double matrix:
%             {a, b} will be become [value(a), value(b)],
%             but {a; b} will remain as {value(a); value(b)}
%
% RETURNS:
%   v       - Value with no AD status.
%
% SEE ALSO:
%   double2ADI, double2GenericAD, ADI, GenericAD

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
    if nargin < 2
        flatten = true;
    end
    if isnumeric(x) || islogical(x) || ischar(x)
        xv = x;
    elseif iscell(x)
        sz = size(x);
        xv = applyFunction(@(x) value(x, flatten), x);
        if flatten && sz(1) == 1
            % Cell arrays are converted to matrices
            xv = [xv{:}];
        end
    elseif isstruct(x)
        fn = fieldnames(x);
        for i = 1:numel(fn)
            f = fn{i};
            for j = 1:numel(x)
                x(j).(f) = value(x(j).(f), flatten);
            end
        end
        xv = x;
    else
        try
            % If the object somehow has the property '.val' we use that. Note
            % that this is partially to fix an Octave bug when using function
            % handles on classes.
            xv = x.val;
        catch
            xv = x;
        end
    end
end
