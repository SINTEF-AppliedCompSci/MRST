function c = getVEColors(ind)
% Return colors used for plotting VE models
%
% SYNOPSIS:
%   c = getVEColors(ind)
%
% PARAMETERS:
%   ind - a scalar or vector consisting of one or more of the numbers one
%         to five, or a cell array consisting of a set of strings, that
%         signify the following parts of a VE model:
%            1 - trapped
%            2 - plume
%            3 - residual
%            4 - dissolved
%            5 - brine
%
% RETURNS:
%   c  - a numel(ind)x3 vector specifying one RGB color per entry in 'ind'
%
% EXAMPLE:
%   image(repmat(1:5,2,1)); colormap(getVEColors(1:5))
%
% SEE ALSO:
%   `getInventoryColors`

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

red = [1 0 0]; blue = [0 0 1]; white = [1 1 1]; black = [0 0 0];
mix = @(a,b,c) (a*c + b*(100-c))/100;

if ischar(ind)
   ind = {ind};
elseif isnumeric(ind)
   ind = num2cell(ind);
end
c = zeros(numel(ind),3);
for i=1:numel(ind)
   switch ind{i}
      case {1,'trapped'}
         c(i,:) = mix(mix(red,white,40),black,80);
      case {2,'plume'}
         c(i,:) = mix(red,white,40);
      case {3,'residual'}
         c(i,:) = mix(mix(red,blue,50),white,20);
      case {4,'dissolved'}
         c(i,:) = mix(blue,white,20);
      case {5,'brine'}
         c(i,:) = mix(blue,white,30);
   end
end 