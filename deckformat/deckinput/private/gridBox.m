function varargout = gridBox(varargin)
%Input box accessor function.
%
% SYNOPSIS:
%       gridBox(b)
%   b = gridBox
%
% PARAMETERS:
%   b - Input box limits.  Must consist of exactly six (positive) integers
%       that define, respectively, lower and upper limits of the I, J, and
%       K axes of an ECLIPSE/FrontSim simulation model.
%
% RETURNS:
%   b - Current input box limits.

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


   persistent BOX

   if nargin == 1,
      b = varargin{1};

      % Sanity check input box limits.
      %
      assert (isnumeric(b), 'Grid box must be numeric.');
      assert (numel(b) == 6, 'Grid box must consist of six items.');
      assert (all(mod(b, 1) == 0), 'Grid box must be integers only.');
      assert (all(b > 0), 'Grid box limits must be positive.');
      assert (all(diff(reshape(b, [2, 3])) >= 0), ...
              'Grid box limits must specify non-empty index range.');

      BOX = reshape(b, 1, []);
   end

   if nargout == 1,
      assert (numel(BOX) == 6, ...
              'Cannot query input box until limits have been defined.');

      varargout{1} = BOX;
   end
end
