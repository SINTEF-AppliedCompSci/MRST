function [yi, der] = interpTableMEX(X, Y, xi, varargin)
%Interpolate a one-dimensional table with MEX acceleration
%
% SYNOPSIS:
%   yi = interpTable(X, Y, xi)
%   yi = interpTable(X, Y, xi, 'pn1', pv1, ...)
%
% PARAMETERS:
%   X       - Nodes at which underlying function y=y(x) is sampled.
%
%   Y       - Values of the underlying function y=y(x) at the nodes, `X`.
%
%   xi      - Evaluation points for new, interpolated, function values.
%
%   method  - (OPTIONAL) Choose the interpolation method used:
%             1 corresponds to a binary search (optimal for a single lookup
%             in ordered data, default, generally fastest if points are
%             given on an uneven grid)
%             2 will use a binning algorithm which can in very particular
%             cases be faster when a large number of lookups is required.
%             3 will assume that X is completely uniform in grid size. Only
%             the first two entries are used. May give wrong results when
%             assumption is violated - no verification is performed.
% RETURNS:
%   yi - Approximate (interpolated/extrapolated) values of the function
%        y=y(x) at the points `xi`.
%
% SEE ALSO:
%   `interpTable`.
%
% NOTE:
%   This routine has less functionality and validation of inputs than
%   interpTable, but can be significantly faster.

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
    
    if isnumeric(xi)
        if nargout > 1
            [yi, der] = mexInterp1(X, Y, xi);
        else
            yi = mexInterp1(X, Y, xi);
        end
    else
        yi = xi;
        [yi.val, der] = mexInterp1(X, Y, xi.val, varargin{:});
        if any(der)
            yi.jac = yi.lMultDiag(der, yi.jac);
        else
            yi.jac = cellfun(@(x) 0*x, yi.jac, 'UniformOutput', false);
        end
    end
end
