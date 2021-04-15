classdef InteractivePoint < handle
    properties
        Callback
        Parent
        moveMarker
        moveCallBack
    end
    properties (Dependent)
        Position
    end
    methods
        function p = InteractivePoint(varargin)
            opt = struct('Callback',    @(src, event)disp('Hello'), ...
                'moveCallback',@(src,event)disp('move'), ...
                'XData', 0, 'YData', 0, ...
                'Parent',     gca, ...
                'Marker',     's', ...
                'MarkerSize',  10, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'm');
            
            [opt, mopt] = merge_options(opt, varargin{:});
            mopt0 = {'Marker',     's', ...
                'MarkerSize',  10, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'm'};
            p.moveMarker = line('Visible', 'off', 'LineStyle', 'none', ...
                'XData', opt.XData, 'YData', opt.YData, 'ButtonDownFcn', ...
                @p.movePoint, 'Parent', opt.Parent, mopt0{:}, mopt{:});
            
            p.Parent    = p.moveMarker.Parent;
            p.Callback  = opt.Callback;
            p.moveCallBack = opt.moveCallback;
            set(p.moveMarker, 'Visible', 'on');
        end
        
        function set.Position(p, val)
            set(p.moveMarker, {'XData', 'YData'}, {val(1), val(2)});
        end
        
        function val = get.Position(p)
            val = [p.moveMarker.XData, p.moveMarker.YData];
        end
        
        function movePoint(p ,src, event)
            if event.Button == 1
                ax = p.Parent;
                freeze(ax);
                fg = ax.Parent;
                cp = ax.CurrentPoint(1,1:2);
                lp = [src.XData, src.YData];
                motionFcn = fg.WindowButtonMotionFcn;
                upFcn     = fg.WindowButtonUpFcn;
                offset = lp-cp;
                fg.WindowButtonMotionFcn = {@dragPoint, src, offset};
                fg.WindowButtonUpFcn     = {@clearMotionFcn, motionFcn, upFcn};
            end
            function dragPoint(src, event, pnt, offset)
                ax.CurrentPoint;
                newp = ax.CurrentPoint(1, 1:2)+offset;
                [xl, yl] = deal(ax.XLim, ax.YLim);
                set(pnt, {'XData', 'YData'}, {max(xl(1), min(newp(1), xl(2))), max(yl(1), min(newp(2), yl(2)))});
                p.moveCallBack(src, event);
            end
            
            function clearMotionFcn(src, event, motionFcn, upFcn)
                src.WindowButtonMotionFcn = motionFcn;
                src.WindowButtonUpFcn     = upFcn; 
                p.Callback(src, event);
                %r.reDraw(pnt.Tag);
                %defreeze(ax);
            end
        end
        
        
        
    end
end


function delete(l)
delete(l.markers);
delete(l.trajectory);
if isvalid(l.Parent)
    l.Parent.UIContextMenu = gobjects(0);
    l.Parent.ButtonDownFcn = '';
    
    if isvalid(l.Parent.Parent)
        l.Parent.Parent.WindowButtonMotionFcn = '';
        l.Parent.Parent.WindowButtonUpFcn = '';
        l.Parent.Parent.Pointer = 'arrow';
    end
end
end


function freeze(ax)
set(ax, {'XLimMode', 'YLimMode', 'ZLimMode'}, {'manual', 'manual', 'manual'});
end

function defreeze(ax)
set(ax, {'XLimMode', 'YLimMode', 'ZLimMode'}, {'auto', 'auto', 'auto'});
end




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

