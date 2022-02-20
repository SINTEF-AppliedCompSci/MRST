classdef InteractiveRectangle < handle
    properties
        Callback
        Parent
        moveMarker
        dragMarkers
        Line
        XLimInner
        YLimInner
    end
    properties (Dependent)
        Position
    end
    methods
        function r = InteractiveRectangle(varargin)
            opt = struct('Position', [0 0 1 1], ...
                         'XLimInner',            [], ...
                         'YLimInner',            [], ...
                         'Callback',  @(src, event)'Hello');
            [opt, other] = merge_options(opt, varargin{:});
            
            cc      = opt.Position;
            cc(3:4) = cc(1:2) + cc(3:4);
            r.Line = line('XData', cc([1 3 3 1 1]'), ...
                          'YData', cc([2 2 4 4 2]'), ...
                          other{:}, 'PickableParts', 'none');
            %r.Rectangle = rectangle(other{:});
            r.Parent    = r.Line.Parent;
            r.Callback  = opt.Callback;
            
            [r.XLimInner, r.YLimInner] = deal(opt.XLimInner, opt.YLimInner);
            
            % interactive markers
            mopt = {'Marker',     'd', ...
                'MarkerSize',  10, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'w'};
            tags = {'00', '01', '11', '10'};
            pnts = cc;
            %pnts(3:4) = pnts(1:2) + pnts(3:4);
            pos  = {[1 2], [3 2], [3 4], [1 4]};
            for k = 1:4
                dragMarkers(k) = line('Visible', 'off', 'LineStyle', 'none', ...
                    'XData', pnts(pos{k}(1)), ...
                    'YData', pnts(pos{k}(2)), ...
                    'ButtonDownFcn', @r.movePoint, ...
                    'Parent', r.Parent, ...
                    'Tag', tags{k}, mopt{:});
            end
            % move marker
            mopt = {'Marker',     's', ...
                'MarkerSize',  10, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'w'};
            moveMarker = line('Visible', 'off', 'LineStyle', 'none', ...
                    'XData', (pnts(1)+pnts(3))/2, ...
                    'YData',  pnts(2), ...
                    'ButtonDownFcn', @r.movePoint, ...
                    'Parent', r.Parent, mopt{:});
            
            r.dragMarkers = dragMarkers;
            r.moveMarker  = moveMarker;
            set(r.dragMarkers, 'Visible', 'on');
            set(r.moveMarker, 'Visible', 'on');
        end
        
        function set.Position(r, val)
            [px,py] = deal([val(1), val(1)+val(3)], [val(2), val(2)+val(4)]);
            set(r.Line, {'XData', 'YData'}, {px([1 2 2 1 1]'), py([1 1 2 2 1]')});
            r.reDraw();
        end
        
        function val = get.Position(r)
            [x,y] = deal(r.Line.XData, r.Line.YData);
            val = [x(1), y(1), x(2)-x(1), y(3)-y(1)];
        end
        function movePoint(r ,src, event)
            if event.Button == 1
                ax = r.Parent;
                freeze(ax);
                fg = ax.Parent;
                cp = ax.CurrentPoint(1,1:2);
                lp = [src.XData, src.YData];
                motionFcn = fg.WindowButtonMotionFcn;
                offset = lp-cp;
                fg.WindowButtonMotionFcn = {@dragPoint, src, offset};
                fg.WindowButtonUpFcn     = {@clearMotionFcn, motionFcn};
            end
            function dragPoint(src, event, pnt, offset)
                ax.CurrentPoint;
                newp = ax.CurrentPoint(1, 1:2)+offset;
                
                %[x1, y1] = deal(max(xl(1), min(newp(1), xl(2))), max(yl(1), min(newp(2), yl(2))));
                [xl, yl] = deal(ax.XLim, ax.YLim);
                set(pnt, {'XData', 'YData'}, {max(xl(1), min(newp(1), xl(2))), max(yl(1), min(newp(2), yl(2)))});
                r.reDraw(pnt.Tag);
            end
            
            function clearMotionFcn(src, event, motionFcn)
                src.WindowButtonMotionFcn = motionFcn;
                r.Callback(src, event)
                %r.reDraw(pnt.Tag);
                defreeze(ax);
            end
        end
        
        function reDraw(r, tag)
            pos = r.Position;
            pos(3:4) = pos(1:2) + pos(3:4);
            if nargin == 2
                switch tag
                    case '00'
                        pos(1:2)   = [r.dragMarkers(1).XData, r.dragMarkers(1).YData];
                    case '01'
                        pos([3 2]) = [r.dragMarkers(2).XData, r.dragMarkers(2).YData];
                    case '11'
                        pos(3:4)   = [r.dragMarkers(3).XData, r.dragMarkers(3).YData];
                    case '10'
                        pos([1 4]) = [r.dragMarkers(4).XData, r.dragMarkers(4).YData];
                    otherwise
                        xshift = r.moveMarker.XData - (pos(1)+pos(3))/2;
                        yshift = r.moveMarker.YData - pos(2);
                        pos([1 3]) =  pos([1 3]) + xshift;
                        pos([2 4]) =  pos([2 4]) + yshift;
                end
            end
            if ~isempty(r.XLimInner)
                if pos(1) > r.XLimInner(1) || pos(3) < r.XLimInner(2)
                    pos = r.Position;
                    pos(3:4) = pos(1:2) + pos(3:4);
                end
            end
            if ~isempty(r.YLimInner)
                if pos(2) > r.YLimInner(1) || pos(4) < r.YLimInner(2)
                    pos = r.Position;
                    pos(3:4) = pos(1:2) + pos(3:4);
                end
            end
            
            set(r.Line, {'XData', 'YData'}, {pos([1 3 3 1 1]'), pos([2 2 4 4 2]')});
            xx = pos([1 3 3 1]);
            yy = pos([2 2 4 4]);
            for k = 1:4
                set(r.dragMarkers(k), {'XData', 'YData'}, {xx(k), yy(k)});
            end
            r.moveMarker.XData = sum(xx(1:2))/2;
            r.moveMarker.YData = yy(1);
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
        
