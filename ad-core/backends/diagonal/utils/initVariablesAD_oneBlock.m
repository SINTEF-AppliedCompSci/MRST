function varargout = initVariablesAD_oneBlock(varargin)
% Initialize a set of automatic differentiation variables (single block)

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

   offsets = cumsum([1, numvals])';
   ncols = sum(numvals);
   for i = 1 : n
      nrows  = numvals(i);
  
      I = 1:nrows;
      J = (1:nrows) + offsets(i) - 1;
      jac = {sparse(I, J, 1, nrows, ncols)};

      varargout{i} = GenericAD(varargin{i}, jac);
      varargout{i}.numVars = numvals';
   end
end
