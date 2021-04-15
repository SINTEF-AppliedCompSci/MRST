function data = readStructField(s, name)
%Read a field from a struct.
%
% SYNOPSIS:
%   data = readStructField(rock, 'rock.perm:3')
%
% DESCRIPTION:
%   The inverse function of getStructFields. Acts in much the same manner
%   as s.(name), but has some additional logic to handle multi-column
%   datasets.
%
% REQUIRED PARAMETERS:
%   s - struct or numerical vector used to produce the name input by
%       getStructFields.
%
%   name - Name of field requested as output. Must be a valid field as
%          defined by getStructFields.
% RETURNS:
%   data - A Nc x 1 vector or Nc x 3 matrix suitable for plotCellData.
%
%
% SEE ALSO:
%   `getStructFields`, `datasetSelector`

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

    f = regexp(name, '\.', 'split');
    jind = 1;
    for i = 2:numel(f)
        fname = f{i};
        subind = regexp(f{i}, ':', 'split');
        if numel(subind) > 1
            fname = subind{1};
            jind = str2double(subind{2});
        else
            jind = 1;
        end
        s = s.(fname);
    end
    if size(s, 2) == 3 && numel(subind) == 1 && isnumeric(s)
        jind = 1:3;
        s = bsxfun(@rdivide, s, sum(s,2));
    end
    if iscell(s)
        data = s{jind};
    else
        data = s(:, jind);
    end
end
