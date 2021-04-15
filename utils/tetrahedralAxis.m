function [mapx, mapy, mapz] = tetrahedralAxis(varargin)
%Create a ternary axis and mappings to ternary space
%
% SYNOPSIS:
%       [mapx, mapy] = ternaryAxis()
%
% OPTIONAL PARAMETERS:
%   'tick'  - Tick points to be displayed on the axis (vector from 0 to 1)
%
%   'names' - Names of the x, y, z coordinates (to be plotted on axis)
%
%   'isox','ixoy','isoz' - Booleans indicating if isolines are to be drawn
%                          for x, y and z axes. Isolines are drawn at the
%                          same positions as the ticks.
%
% RETURNS:
%   mapx  - Mapping on the form f(x, y, z) -> X where X is the new
%           coordinates inside the ternary diagram.
%
%   mapy  - Mapping on the form g(x, y, z) -> Y where Y is the new
%           coordinates inside the ternary diagram.
%
% EXAMPLE:
% figure; 
% [mapx, mapy] = ternaryAxis();
% x = 0.7*(0:0.1:1); 
% y = 0.3*(0:0.1:1).^2;
% z = 1 - x - y; 
% plot(mapx(x, y, z), mapy(x, y, z), 'k', 'linewidth', 2)
%
% SEE ALSO:
%   `ternaryAxis`

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


    opt = struct('names', {{'X', 'Y', 'Z', 'L'}} ...
                 );
    opt = merge_options(opt, varargin{:});
    hold on
    
    mapx = @(x, y, z, l) (1/2)*(z + 1 - y)./(x + y+ z + l) - 1/4;
    mapy = @(x, y, z, l) sqrt(3)*(x/2 + l/6)./(x + y+ z + l) - sqrt(3)*(1/8 + 1/24);
    mapz = @(x, y, z, l) (sqrt(6)/3).*l./(x + y+ z + l) - (sqrt(6)/12);
    
    
    x0 = [mapx(1, 0, 0, 0), mapy(1, 0, 0, 0), mapz(1, 0, 0, 0)];
    y0 = [mapx(0, 1, 0, 0), mapy(0, 1, 0, 0), mapz(0, 1, 0, 0)];
    z0 = [mapx(0, 0, 1, 0), mapy(0, 0, 1, 0), mapz(0, 0, 1, 0)];
    l0 = [mapx(0, 0, 0, 1), mapy(0, 0, 0, 1), mapz(0, 0, 0, 1)];
    
    plot3([x0(1); y0(1); z0(1); x0(1)], ...
          [x0(2); y0(2); z0(2); x0(2)], ...
          [x0(3); y0(3); z0(3); x0(3)], 'k')

    plot3([x0(1); y0(1); l0(1); x0(1)], ...
          [x0(2); y0(2); l0(2); x0(2)], ...
          [x0(3); y0(3); l0(3); x0(3)], 'k')

    plot3([x0(1); z0(1); l0(1); x0(1)], ...
          [x0(2); z0(2); l0(2); x0(2)], ...
          [x0(3); z0(3); l0(3); x0(3)], 'k')
 
    plot3([y0(1); l0(1); z0(1); y0(1)], ...
          [y0(2); l0(2); z0(2); y0(2)], ...
          [y0(3); l0(3); z0(3); y0(3)], 'k')
      
    if ~isempty(opt.names)
        text(x0(1), x0(2), x0(3), [opt.names{1},'=1'], 'verticalalignment', 'top', 'horizontalalignment', 'center')
        text(y0(1), y0(2), y0(3), [opt.names{2},'=1'], 'verticalalignment', 'top', 'horizontalalignment', 'center')
        text(z0(1), z0(2), z0(3), [opt.names{3},'=1'], 'verticalalignment', 'top', 'horizontalalignment', 'center')
        text(l0(1), l0(2), l0(3), [opt.names{4},'=1'], 'verticalalignment', 'bottom', 'horizontalalignment', 'center')
    end
    

    axis equal tight off
end
