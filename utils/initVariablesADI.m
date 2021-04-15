function varargout = initVariablesADI(varargin)
% Initialize a set of automatic differentiation variables
%
% SYNOPSIS:
%  a            = initVariablesADI(a);
%  [a, b, c, d] = initVariablesADI(a, b, c, d);
%
% PARAMETERS:
%   varargin - Any number of variables in either column vector format or as
%              scalars. These variables will be instantiate as ADI objects
%              containing both a .val field and a .jac jacobian. These
%              variables will start with identity jacobians with regards to
%              themselves and zero jacobians with regards to the other
%              variables (implicitly defined by the ordering of input and
%              output).
%
%              These variables can then be used to create more complex
%              expressions, resulting in automatic compuation of the first
%              order derivatives leading to easy implementation of
%              Newton-like nonlinear solvers.
%
% EXAMPLE:
%        x = 1;
%        y = 5;
%        [x, y] = initVariablesADI(x, y)
%
%        This gives x.jac ->  {[1]  [0]} and y.jac ->  {[0]  [1]}.
%
%        If we compute z = x.*y.^2 we get
%
%        z.val = 25 (as is expected),
%        z.jac{1} = d(x*y^2)/dx = y^2 = 5^2 = 25
%        z.jac{2} = d(x*y^2)/dy = 2*x*y = 2*1*5 = 10;
%
%        Note that as this is meant for vector operations, the
%        element-wise operations should be used (.* instead of *) even when
%        dealing with scalars.
%
% RETURNS:
%   varargout - The same variables as inputted, as ADI objects.
%
% SEE ALSO:
%   `ADI`

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


   assert (nargin == nargout, ...
          ['Number of output variables must equal the number ', ...
           'of input variables.']);

   numvals   = reshape(cellfun('prodofsize', varargin), 1, []);
   n         = nargin;
   varargout = cell([1, n]);

   for i = 1 : n,
      j = [ 1 : (i - 1), (i + 1) : n ];

      nrows  = numvals(i);
      jac    = cell([1, n]);
      jac(j) = arrayfun(@(ncols) sparse(nrows, ncols), ...
                        numvals(j), 'UniformOutput', false);

      jac{i} = sparse(1 : nrows, 1 : nrows, 1, nrows, nrows);

      varargout{i} = ADI(varargin{i}, jac);
   end
end
