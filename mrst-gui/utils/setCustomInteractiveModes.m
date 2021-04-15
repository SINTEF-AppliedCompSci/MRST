function s = setCustomInteractiveModes(obj)
% set custom props for zoom/pan/rotate in axes ax

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

% apply to given axis (and figure parent

% enable panel-clicks while in interacitve mode
s = struct('rotate3d', [], 'zoom', [], 'pan', []);

if strcmp(obj.Type, 'axes')
    ax  = obj;
    fig = ax.Parent;
elseif strcmp(obj.Type, 'figure')
    fig = obj;
    c   = fig.Children;
    ii  = strcmp(get(c, 'Type'), 'axes');
    ax  = c(ii);
end

if isempty(ax)
    warning('Setting interactive modes for figure without axes will only have partial effect');
end

s.zoom = zoom(fig);
s.zoom.ButtonDownFilter = @filterfunc;
s.pan = pan(fig);
s.pan.ButtonDownFilter = @filterfunc;
s.rotate3d = rotate3d(fig);
s.rotate3d.ButtonDownFilter = @filterfunc;
% set mode context menus
setMenus = true;
if isprop(ax(1), 'ContextMenu')
    cmenu = 'ContextMenu';
elseif isprop(ax(1), 'UIContextMenu')
    cmenu = 'UIContextMenu';
else
    warning('Unable to set axes/zoom/pan/ratate context menus');
    setMenus = false;
end

if setMenus
    s.zoom.(cmenu)     = getModeContextMenu(s, ax, 'zLabel', 'Zoom off');
    s.pan.(cmenu)      = getModeContextMenu(s, ax, 'pLabel', 'Pan off');
    s.rotate3d.(cmenu) = getModeContextMenu(s, ax, 'rLabel', 'Rotate off');
    % set axes context menus/style
    for k = 1:numel(ax)
        ax(k).(cmenu) = getModeContextMenu(s, ax);
        if ~verLessThan('matlab', '9.1')
            setAxes3DPanAndZoomStyle(s.zoom, ax(k), 'camera');
        end
    end
end
end

%--------------------------------------------------------------------------
function m = getModeContextMenu(s, ax, varargin)
opt = struct('rLabel', 'Rotate', ...
    'zLabel', 'Zoom', ...
    'pLabel', 'Pan');
opt = merge_options(opt, varargin{:});
m = uicontextmenu('Parent', ax(1).Parent);
uimenu('Parent', m, 'Label', opt.rLabel', 'Callback', @toggleMode, 'Tag', 'rotate3d');
uimenu('Parent', m, 'Label', opt.zLabel, 'Callback', @toggleMode, 'Tag', 'zoom');
uimenu('Parent', m, 'Label', opt.pLabel, 'Callback', @toggleMode, 'Tag', 'pan');
uimenu('Parent', m, 'Label', 'Reset', 'Callback', @resetZoom, 'Tag', 'reset')
    function toggleMode(src, ~)
        if strcmp(s.(src.Tag).Enable, 'off')
            s.(src.Tag).Enable = 'on';
        else
            s.(src.Tag).Enable = 'off';
        end
    end
    function resetZoom(~, ~)
        if ~isempty(ax)
            for k = 1:numel(ax)
                zoom(ax(k), 'reset')
            end
        end
    end
end

%--------------------------------------------------------------------------
function bol = filterfunc(src, ~)
%bol = isa(src,'matlab.ui.container.Panel') || isa(src,'matlab.ui.Figure')
bol = ~any(strcmp(src.Type, {'axes', 'patch', 'line'}));
end


