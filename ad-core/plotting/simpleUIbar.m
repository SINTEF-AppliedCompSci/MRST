function simpleUIbar(parent, data, start, height, txt, varargin)
%Undocumented Utility Function

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

    h = height/2;
    set(0, 'CurrentFigure', parent);
    axes('position', [0.1, start + h, .8, h])
    plotProgress(data, varargin{:})
    
    uicontrol(parent, 'units', 'Normalized', ...
        'Style', 'text', 'string', txt,...
        'position', [0.1, start, .8, h*0.9])
end

function plotProgress(data, varargin)
    rectangle('Position', [0, 0.01, data, 0.99], 'FaceColor', 'r', varargin{:})
    hold on
    rectangle('Position', [0, 0.01, 1, 0.99], 'FaceColor', 'none')
    ylim([0, 1]);
    xlim([0, 1]);
    axis off
    txt = sprintf('%2.1f%%', 100*data);
    text(0.5, 0.5, txt)
end
