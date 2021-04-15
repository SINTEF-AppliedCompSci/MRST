function [fx, cx] = splitFaceCellValue(operators, flag, x, sz)
% Split multi-valued function into cell and face values
%
% SYNOPSIS:
%   [fx, cx] = splitFaceCellValue(operators, flag, x, sz)
%
% DESCRIPTION:
%   Taking a set of values, this function returns cell and face-upwinded
%   values based on specified flag and dimensions. Normally, this is simply
%   applying a pre-existing upwind operator to get the upstream weighted
%   values for transported quantities. For special functions that arise in
%   some workflows, it can take e.g. a set of (half)face values plus cell
%   values and divide them up in a reasonable manner.
%
% PARAMETERS:
%   operators - Operators struct. See `setupOperatorsTPFA`.
%   flag      - Upstream flag to be used to upwind values. See `faceUpstr`.
%   x         - Vector of values to be treated. Can be either one value per
%               cell in the domain, one value per face followed by one
%               value per cell, or one value per half-face, followed by the
%               cell values.
%   sz        - Vector of length 2. First entry corresponds to the number
%               of faces and the second is the total number of cells.
%
% SEE ALSO:
%   `faceUpstr`, `setupOperatorsTPFA`

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

    nf = sz(1);
    nc = sz(2);
    n = numelValue(x);

    switch n
        case nc
            % Cell-wise values only, use upstream weighting
            fx = operators.faceUpstr(flag, x);
            cx = x;
        case nf + nc
            % Face values first, then cell values
            fx = x(1:nf);
            cx = x((nf+1):end);
        case 2*nf + nc
            % Half face values
            subs = (1:nf)' + ~flag.*nf;
            fx = x(subs);
            cx = x((2*nf+1):end);
        case nf
            % Only face values
            error('Not implemented yet');
        case 2*nf
            % Only half-face values
            error('Not implemented yet');
        otherwise
            error('Did not find expected dimension of input');
    end
end