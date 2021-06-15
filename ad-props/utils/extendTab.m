function T = extendTab(T, pm)
   if nargin < 2, pm = 1; end

   if ~ iscell(T),
      T = extend(T, pm);
   else
      T = cellfun(@(t) extend(t, pm), T, 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function T = extend(T, pm)
   t1 = T( 1 ,:); t1(1) = t1(1) - pm;
   te = T(end,:); te(1) = te(1) + pm;

   T = [t1; T; te];
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
