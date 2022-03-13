function c = getInventoryColors(ind)
% Return colors used for plotting CO2 inventories (green to orange color)
%
% SYNOPSIS:
%   c = getInventoryColors(ind)
%
% PARAMETERS:
%   ind - a scalar or vector consisting of one or more of the numbers one
%         to seven. These numbers signify the following parts of a CO2
%         inventory:
%            1 - dissolved
%            2 - residual (traps)
%            3 - residual
%            4 - residual (plume)
%            5 - movable (traps)
%            6 - movable (plume)
%            7 - leaked
%
% RETURNS:
%   c  - a numel(ind)x3 vector specifying one RGB color per entry in 'ind'
%
% EXAMPLE:
%   image(repmat(1:7,2,1)); colormap(getInventoryColors(1:7))
%
% SEE ALSO:
%   `interactiveTrapping`, `getVEColors`

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
red    = [1 0 0];
green  = [0 1 0];
yellow = [1 1 0];
orange = [1 .65 0];
white  = [1 1 1]; 
black  = [0 0 0];
mix = @(a,b,c) (a*c + b*(100-c))/100;

c = zeros(numel(ind),3);
for i=1:numel(ind)
   switch ind(i)
      case 1
         c(i,:) = mix(green,black,50);
      case 2
         c(i,:) = mix(green,black,75);
      case 3
         c(i,:) = green;
      case 4
         c(i,:) = mix(green,white,50);
      case 5
         c(i,:) = mix(yellow,white,80);
      case 6
         c(i,:) = mix(orange,white,80);
      case 7
         c(i,:) = mix(orange,red,30);
   end
end 
