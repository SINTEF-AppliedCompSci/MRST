function varargout = defaultBox(varargin)
%Default input box accessor function.
%
% SYNOPSIS:
%       defaultBox(d)
%   b = defaultBox
%
% PARAMETERS:
%   d - Input box limits.  Must consist of exactly three (positive)
%       integers that define the default upper input box limits of the I,
%       J, and K axes of an ECLIPSE/FrontSim simulation model.  These
%       limits are assumed to correspond to the entire model and should
%       equal the data input in the 'DIMENS' keyword (i.e., the 'cartDims'
%       field of the MRST grid_structure.)
%
% RETURNS:
%   b - Default input box limits defined by [1, d(1), 1, d(2), 1, d(3)].
%
% SEE ALSO:
%   `grid_structure`, `private/gridBox`.

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


   persistent DIMENS

   if nargin == 1,
      d = varargin{1};

      % Sanity check input box limits.
      %
      assert (isnumeric(d), 'Grid box must be numeric.');
      assert (numel(d) == 3, 'Grid box must consist of three items.');
      assert (all(mod(d, 1) == 0), 'Grid box must be integers only.');
      assert (all(d > 0), 'Default grid box limits must be positive.');

      DIMENS = reshape(d, 1, []);
   end

   if nargout == 1,
      assert (numel(DIMENS) == 3, ...
             ['Cannot query default input box until upper limits ', ...
              'have been set.']);

      varargout{1} = reshape([ones([1, 3]); DIMENS], 1, []);
   end
end
