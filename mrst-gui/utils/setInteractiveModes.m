function s = setInteractiveModes(ax)
% set props for zoom/pan/rotate in axes ax
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
bol = ~any(strcmp(src.Type, {'axes', 'patch'}));
end