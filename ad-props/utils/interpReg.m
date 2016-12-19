function [yi, dyidxi] = interpReg(T, xi, reginx)
   
   if isempty(xi)
      yi     = [];
      dyidxi = [];
      return;
   end
   
   nreg = numel(reginx);
   
   if nreg > 0 && ischar(reginx{1}) && strcmp(reginx{1}, ':'),
      % Special case denoting entire domain in single region.

      yi = fastInterpTable(T{1}(:,1), T{1}(:,2), xi);
   
   elseif nreg > 0
      
      yi = zeros(size(xi));

      for k = 1:nreg,
         if ~isempty(reginx{k})
            yi(reginx{k}) = fastInterpTable(T{k}(:,1), T{k}(:,2), ...
                                        xi(reginx{k}));
         end
      end
   
   end
   
   if nargout > 1,

      dyidxi = dinterpReg(T, xi, reginx);

   end
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
