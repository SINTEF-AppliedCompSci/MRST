function x = geomspace(a, b, n, x0)
%Geometrically spaced vector.
%
% SYNOPSIS:
%   x = geomspace(x1, x2, n, L0)
%
% PARAMETERS:
%   x1,x2  - Lower and upper bounds on resulting vector.
%
%   n      - Number of points to generate between `x1` and `x2`, inclusive.
%            Must be at least two (`n` >= 2).
%
%   L0     - Length of first (and shortest) sub-interval.
%            Must not exceed `(x2 - x1) / (n - 1)`.
%
% RETURNS:
%   x - n-point row vector from x1 to x2 constructed such that the length
%       of consecutive sub-intervals (i.e., DIFF(x)) *increases* by a
%       constant (geometric) factor.
%
% NOTE:
%   This function is based on solving a polynomial equation of degree n-1
%   (using the `roots` function).  This is an unstable process and,
%   consequently, `n` should not be too big (usually no greater than 100).
%
%   Due to round-off errors, x(end) may differ slightly from x2.
%
%   These facts mean that geomspace is not a general utility and its use
%   should be carefully considered in every instance.
%
% SEE ALSO:
%   `linspace`, `logspace`, `roots`.

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


   assert (n >= 2);

   f = (b - a) / x0;
   if ~(f > n - 1),
      error(msgid('GeomSequence:DoesNotExist'), ...
            'No geometric sequence with these parameters exists.')
   end

   c        = zeros([n, 1]);
   c(    1) = 1;
   c(end-1) = -f;
   c(end  ) = f - 1;

   q = roots(c);
   q = q(~(abs(imag(q)) > 0) & (real(q) > 1));   assert (numel(q) == 1);

   x = cumsum([a, x0 .* q.^(0 : n-2)]);
end
