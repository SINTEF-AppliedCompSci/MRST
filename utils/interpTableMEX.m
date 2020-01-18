function yi = interpTableMEX(X, Y, xi, varargin)
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
%
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
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
        yi = mexInterp1(X, Y, xi);
    else
        yi = xi;
        [yi.val, der] = mexInterp1(X, Y, xi.val);
        yi.jac = yi.lMultDiag(der, yi.jac);
    end
end
