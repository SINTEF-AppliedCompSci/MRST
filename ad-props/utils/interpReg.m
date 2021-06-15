function [yi, dyidxi] = interpReg(T, xi, reginx)
%Interpolate data (with region support)
%
% SYNOPSIS:
%   [yi, dyidxi] = interpReg(T, xi, reginx)
%
% REQUIRED PARAMETERS:
%   T      - Interpolation table. Cell array of N tables for interpolation
%            where N is the number of regions. Table should consist of a
%            function f(x) in the format [x, f(x)] where x are points and 
%            f(x) is the function values, both in the format of column
%            vectors.
%
%   xi     - Query points where the function should be interpolated.
%
%   reginx - Cell array of length N, where each entry corresponds to a list
%            of the cells present in a certain region. See getRegMap for
%            details.
%
% RETURNS:
%   yi     - Function evaluated at given points, accounting for regions.
%
%   dyidxi - Slope at xi (only computed if requested).
%

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

   if isempty(xi)
      yi     = [];
      dyidxi = [];
      return;
   end
   nreg = numel(reginx);
   
   if nreg > 0 
      if ischar(reginx{1}) && strcmp(reginx{1}, ':'),
          % Special case denoting entire domain in single region.
          yi = fastInterpTable(T{1}(:,1), T{1}(:,2), xi);
      else
          % Iterate over region indices and assign values
          yi = zeros(size(xi));
          for k = 1:nreg,
             subs = reginx{k};
             if ~isempty(subs)
                yi(subs) = fastInterpTable(T{k}(:,1), T{k}(:,2), xi(subs));
             end
          end
      end
   end
   % Derivatives was requested, compute them
   if nargout > 1,
      dyidxi = dinterpReg(T, xi, reginx);
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
