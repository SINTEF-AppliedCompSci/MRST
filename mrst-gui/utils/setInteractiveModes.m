function s = setInteractiveModes(ax)
% set props for zoom/pan/rotate in axes ax

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

% enable panel-clicks while in interacitve mode
fig = ax.Parent;

s.zoom = zoom(fig);
s.zoom.ButtonDownFilter = @filterfunc;

s.pan = pan(fig);
s.pan.ButtonDownFilter = @filterfunc;

s.rotate3d = rotate3d(fig);
s.rotate3d.ButtonDownFilter = @filterfunc;

if ~verLessThan('matlab', '9.1')
    setAxes3DPanAndZoomStyle(s.zoom,ax,'camera');
end
end

function bol = filterfunc(src, ~)
%bol = isa(src,'matlab.ui.container.Panel') || isa(src,'matlab.ui.Figure')
bol = ~any(strcmp(src.Type, {'axes', 'patch', 'line'}));
end
