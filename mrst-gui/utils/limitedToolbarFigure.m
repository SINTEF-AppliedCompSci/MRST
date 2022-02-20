function f = limitedToolbarFigure(varargin)
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

f = figure(varargin{:});
f.MenuBar = 'none';
f.ToolBar = 'figure';
f.DockControls = 'off';

cl = allchild(allchild(f));

% list which buttons to keep and not:
% keep  -  tooltipstring
l = {...
0, 'Show Plot Tools and Dock Figure',...
0, 'Hide Plot Tools',...
0, 'Insert Legend',...
0, 'Insert Colorbar',...
0, 'Link Plot',...
0, 'Brush/Select Data',...
0, 'Data Cursor',...
1, 'Rotate 3D',...
1, 'Pan',...
1, 'Zoom Out',...
1, 'Zoom In',...
0, 'Edit Plot',...
0, 'Print Figure',...
0, 'Save Figure',...
0, 'Open File',...
0, 'New Figure'};

l = l(2*find(cell2mat(l(1:2:end))));
set(cl, 'Visible', 'off');
set(cl, 'Separator', 'off');
for k = 1:numel(l)
    c = findall(cl,'ToolTipString', l{k});
    set(c, 'Visible', 'on');
end
end
