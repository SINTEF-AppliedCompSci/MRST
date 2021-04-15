function [mapx, mapy] = ternaryAxis(varargin)
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
%   `tetrahedralAxis`

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


    opt = struct('tick', 0.2:0.2:0.8, ...
                 'tickLabels', true, ...
                 'names', {{'X', 'Y', 'Z'}}, ...
                 'isox', true, ...
                 'isoy', true,...
                 'isoz', true...
                 );
    [opt, extra] = merge_options(opt, varargin{:});
    hold on
    
    mapx = @(x, y, z) (1/2)*(2*y + z)./(x + y+ z);
    mapy = @(x, y, z) (sqrt(3)/2)*z./(x + y+ z);

    if ~isempty(opt.names)
        text(mapx(1, 0, 0), mapy(1, 0, 0), [opt.names{1},'=1'], 'verticalalignment', 'top', 'horizontalalignment', 'center', extra{:})
        text(mapx(0, 1, 0), mapy(0, 1, 0), [opt.names{2},'=1'], 'verticalalignment', 'top', 'horizontalalignment', 'center', extra{:})
        text(mapx(0, 0, 1), mapy(0, 0, 1), [opt.names{3},'=1'], 'verticalalignment', 'bottom', 'horizontalalignment', 'center', extra{:})
    end
    plot([0, 0.5, 1, 0], [0, sqrt(3)/2, 0, 0], 'k')

    axis off
    
    if ~isempty(opt.tick)
        for i = 1:numel(opt.tick)
            ti = opt.tick(i);
            if opt.tickLabels
                txt = num2str(ti);
                text(mapx(1-ti, ti, 0), mapy(1-ti, ti, 0), txt, ...
                    'verticalalignment', 'top', 'horizontalalignment', 'center', extra{:})
                text(mapx(0, 1-ti, ti), mapy(0, 1-ti, ti), txt, ...
                    'verticalalignment', 'bottom', 'horizontalalignment', 'left', extra{:})
                text(mapx(ti, 0.0, 1-ti), mapy(ti, 0.0, 1-ti), txt, ...
                    'verticalalignment', 'bottom', 'horizontalalignment', 'right', extra{:})
            end
            seg = 0:0.01:1;
            if opt.isox
                x = ti;
                y = (1-ti).*seg;
                z = (1-ti).*(1-seg);
                plot(mapx(x, y, z), mapy(x, y, z), 'k');
            end
            
            if opt.isoy
                x = (1-ti).*seg;
                y = ti;
                z = (1-ti).*(1-seg);
                plot(mapx(x, y, z), mapy(x, y, z), 'k');
            end
            
            
            if opt.isoy
                x = (1-ti).*seg;
                z = ti;
                y = (1-ti).*(1-seg);
                plot(mapx(x, y, z), mapy(x, y, z), 'k');
            end
        end
    end
end
