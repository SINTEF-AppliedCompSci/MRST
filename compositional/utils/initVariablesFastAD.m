function varargout = initVariablesFastAD(varargin)
% Initialise FastAD objects used internally in some compositional code
   assert (nargin == nargout, ...
          ['Number of output variables must equal the number ', ...
           'of input variables.']);

   n_el = cellfun(@numel, varargin);
   m = n_el(1);
   assert(all(n_el == m));
   
   n         = nargin;
   varargout = cell([1, n]);

   for i = 1 : n,
      jac = zeros(m, n);
      jac(:, i) = 1;
      varargout{i} = FastAD(varargin{i}, jac);
   end
end

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
