function dyidx = dinterpq1(x, y, xi)
%Compute derivative of piecewise linear interpolant.
%
% SYNOPSIS:
%   dyi = dintrpq1(x, y, xi)
%
% PARAMETERS:
%   x   - Vector of length at least two containing coordinates of
%         underlying interval.  The values must either be monotonically
%         increasing or monotonically decreasing.
%
%   y   - Function values of piecewise linear function at points x.  Must
%         be a vector containing `numel(x)` elements.
%
%   xi  - Points at which to compute derivative values of function `y(x)`.
%
% RETURNS:
%   dyi - Column vector containing derivative values of function y(x).
%         One scalar value for each point in the input vector `xi`.
%
% NOTE:
%   Function dinterpq1 assumes that all input arrays contain only finite
%   values.  That is, we assume that ::
%
%       ~(any(isinf(x)) || any(isinf(y)) || any(isinf(xi)))
%
%   Furthermore dinterpq1 employs linear extrapolation outside the data
%   points (x,y).
%
% EXAMPLES:
%   % Compute derivatives at xi=[-10, 0.1, pi, 100] of a piecewise linear
%   % function whose piecewise derivative values are 1:10.
%   x = 0 : 10; y = cumsum(x);
%   dinterpq1(x, y, [-10, 0.1, pi, 100])
%
% SEE ALSO:
%   `interp1q`.

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


   assert (numel(x) > 1);
   assert (~any(isinf(xi)));
   assert (~any(isinf(x) |  isinf(y)));
   assert ( all(size (x) == size (y)));

   x = reshape(x, [], 1);
   y = reshape(y, [], 1);

   if any(diff(x) < 0), x = flipud(x); y = flipud(y); end

   dyidx = diff(y) ./ diff(x);     % Derivatives in [min(x),max(x)].
   dyidx = dyidx([1, 1:end, end]); % Linear extrapolation outside interval.

   % Extend interval to cover entire real line.  This ensures b~=0 for all
   % elements of the array xi.  HISTC does not distinguish between the
   % cases xi<min(x) and xi>max(x) (bin is 0 in both cases), but it is
   % important that we be able to compute the correct linear extrapolant
   % for the data.
   %
   [b,b] = histc(xi, [-inf; x; inf]);
   dyidx = reshape(dyidx(b), [], 1);    % Extract correct derivatives.
end
