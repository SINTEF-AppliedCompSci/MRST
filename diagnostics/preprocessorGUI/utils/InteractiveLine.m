classdef InteractiveLine < handle
    properties
        Callback
        Parent
        trajectory
        markers
        smoothing = 'spline';
        nSmooth = 10;
        lineOpt   = {};
        markerOpt = {};
        drawMode  = 'afterLast';
    end
     properties (Dependent)
        XData
        YData
        ZData
     end
    
    methods
        function l = InteractiveLine(varargin)
            if mod(nargin, 2) == 1
                l.smoothing = varargin{1};
                varargin    = varargin(2:end);
            end
            opt       = struct('Parent', [], 'Callback', @(src, event)'');
            [opt, other] = merge_options(opt, varargin{:});
            
            pnts = struct('XData', [], 'YData', [], 'ZData', []);
            [pnts, other] = merge_options(pnts, other{:});
            l.Callback = opt.Callback;
            mopt = struct('Marker',     'o', ...
                          'MarkerSize',  10, ...
                          'MarkerEdgeColor', 'k', ...
                          'MarkerFaceColor', 'g');
            [mopt, other] = merge_options(mopt, other{:});
                           
            mopt = [fieldnames(mopt), cellfun(@(x)mopt.(x), fieldnames(mopt), 'UniformOutput', false)];
            
            l.markerOpt = reshape(mopt', 1, []);
            
            % check if input is on short-form
            iss = cellfun(@ischar, other);
            fs  = find(iss, 1, 'first');
            if fs > 1
                if mod(numel(iss)-fs, 1)
                    l.lineOpt = other(fs:end);
                else
                    l.lineOpt = other(fs+1:end);
                end
            else
                l.lineOpt   = other;
            end
            
            if isempty(opt.Parent)
                l.Parent = gca;
            else
                l.Parent = opt.Parent;
            end
            
            l.trajectory  = line('Visible', 'off', 'Marker', 'none', ...
                                 'DeleteFcn', @(src, event)l.delete, other{:}, ...
                                 'Parent', l.Parent, 'PickableParts', 'none');
            %if ~isempty(pnts.XData)
                l.trajectory.XData = pnts.XData;
                l.trajectory.YData = pnts.YData;
            %end
            if ~isempty(pnts.ZData)
                l.trajectory.ZData = pnts.ZData;
            end
            
            mz = 1;
            if isempty(l.trajectory.ZData)
                mz = [];
            end
            for k = 1:numel(l.trajectory.XData)
                markers(k) = line('Visible', 'off',  l.markerOpt{:}, 'LineStyle', 'none', ...
                                    'XData', l.trajectory.XData(k), ...
                                    'YData', l.trajectory.YData(k), ...
                                    'ZData', l.trajectory.ZData(k*mz), ...
                                    'ButtonDownFcn', @l.movePoint, ...
                                    'Parent', l.Parent);
                addPointContextMenu(markers(k), l)
            end
            if numel(l.trajectory.XData) == 0
                %markers = repmat(line, 0,0);
                markers = [];
            end
            l.markers = markers;
            
            l.resetLinePoints();
            
            l.addAxisContextMenu();
            
            %l.trajectory.ButtonDownFcn = @l.lineCallback;
            %l.line.UIcontextmenu
            %l.markers.ButtonDownFcn = @l.markerCallback;
            %l.line.UIcontextmenu
            
            l.trajectory.Visible = 'on';
            set(l.markers, 'Visible', 'on');
        end
        
        function val = get.XData(l)
            if ~isempty(l.markers)
                val = vertcat(l.markers.XData);
            else
                val = [];
            end
        end
        
        function val = get.YData(l)
            if ~isempty(l.markers)
                val = vertcat(l.markers.YData);
            else
                val = [];
            end
        end
        
        function val = get.ZData(l)
            if ~isempty(l.markers)
                val = vertcat(l.markers.ZData);
            else
                val = [];
            end
        end
        
        function insertPoint(l, p, n)
            nd =  numel(l.XData);
            if nargin < 3
                n =nd + 1;
            end
            
            assert(n>=1 && n <= nd+1);
            
            if numel(p) == 2
                p3 = [];
            else
                p3 = p(3);
            end
            
            h = line(l.markerOpt{:}, 'LineStyle', 'none', ...
                                     'XData', p(1), 'YData', p(2), 'ZData', p3, ...
                                     'ButtonDownFcn', @l.movePoint);
            addPointContextMenu(h, l)
            h.MarkerFaceColor = 'r';
            if ~strcmp(l.drawMode, 'along')
                l.markers(min(n, nd)).MarkerFaceColor = 'g';
            end
            l.markers = [l.markers(1:n-1), h, l.markers(n:nd)];
            l.resetLinePoints();
        end
        
        function deletePoint(l, n)
            nd =  numel(l.XData);
            assert(n>=1 && n <= nd);
            if nargin < 2
                n =nd;
            end
            
            delete(l.markers(n));
            l.markers = [l.markers(1:n-1), l.markers(n+1:nd)];
            l.resetLinePoints();
        end
        
        function cleanUpPoints(l)
            ix = cellfun(@isempty, {l.markers.XData});
            if any(ix)
                l.markers = l.markers(~ix);
            end
        end
        
        function movePoint(l ,src, event)
            if event.Button == 1
                ax = l.Parent;
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
                newp = ax.CurrentPoint(1, 1:2)+offset;
                [xl, yl] = deal(ax.XLim, ax.YLim);
                set(pnt, {'XData', 'YData'}, {max(xl(1), min(newp(1), xl(2))), max(yl(1), min(newp(2), yl(2)))}); 
                l.resetLinePoints;
            end
            
             function clearMotionFcn(src, ~, motionFcn)
                src.WindowButtonMotionFcn = motionFcn;
                l.resetLinePoints;
                %defreeze(ax);
            end
        end
        
        
        function resetLinePoints(l)
            sp = 1:numel(l.XData);
            si = (1:1/l.nSmooth:numel(sp));
          
            if numel(sp) <= 1
                p = [l.XData, l.YData, l.ZData];
            else
                smth = l.smoothing;         
                if numel(l.markers) <= 2
                    smth = 'linear';
                end
                p = interp1(sp, [l.XData, l.YData, l.ZData], si, smth);
            end
            if ~isempty(l.ZData)
                p3 = p(:,3);
            else
                p3 = [];
            end
            if ~isempty(p)
                set(l.trajectory, {'XData', 'YData', 'ZData'}, {p(:,1), p(:,2), p3});
                l.Callback();
            else
                %l.delete();
            end
        end
        
        function addAxisContextMenu(l)
            ax = l.Parent;
            m  = uicontextmenu('Parent', ax.Parent);
            ax.UIContextMenu = m;
            top = uimenu('Parent', m, 'Label', 'Draw points');
            
            %m1 = uimenu('Parent', top, 'Label', 'After last',   'Callback', @(src,event)l.updateMode(src, event, 'afterLast'));
            %m2 = uimenu('Parent', top, 'Label', 'Before first', 'Callback', @(src,event)l.updateMode(src, event, 'beforeFirst'));
            m3 = uimenu('Parent', top, 'Label', 'Along line', 'Callback', @(src,event)l.updateMode(src, event, 'along'));
            m4 = uimenu('Parent', top, 'Label', 'off', 'Callback', @(src,event)l.updateMode(src, event, 'off'));
        end
        
        
        function updateMode(l, src, event, mode)
            l.drawMode = mode;
            switch mode
                case 'off'
                    l.Parent.ButtonDownFcn = '';
                    l.Parent.Parent.Pointer = 'arrow';
                    set(l.markers, 'MarkerFaceColor', 'g');
                case 'afterLast'
                    l.Parent.ButtonDownFcn  = @(src, event)l.putPoint(src, event);
                    l.Parent.Parent.Pointer = 'circle';
                    l.markers(end).MarkerFaceColor = 'r';
                case 'beforeFirst'
                    l.Parent.ButtonDownFcn  = @(src, event)l.putPoint(src, event);
                    l.Parent.Parent.Pointer = 'circle';
                    l.markers(1).MarkerFaceColor = 'r';
                case 'along'
                    l.Parent.ButtonDownFcn  = @(src, event)l.putPoint(src, event);
                    l.Parent.Parent.Pointer = 'circle';
                    set(l.markers, 'MarkerFaceColor', 'r');
            end
                    
%                 
%             if ~strcmp(mode, 'off')
%                 l.Parent.ButtonDownFcn  = @(src, event)l.putPoint(src, event);
%                 l.Parent.Parent.Pointer = 'circle'; 
%             else
%                 l.Parent.ButtonDownFcn = '';
%                 l.Parent.Parent.Pointer = 'arrow'; 
%             end
        end
        
        function putPoint(l, src, event)
            %ax = l.Parent;
            %p = ax.CurrentPoint(1,1:2);
            p = event.IntersectionPoint(1:2);
            if event.Button == 1
                switch l.drawMode
                    case 'afterLast'
                        l.insertPoint(p);
                    case 'beforeFirst'
                        l.insertPoint(p, 1);
                    case 'along'
                        ix = getPointInsertIx(l, p);
                        l.insertPoint(p, ix);
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

        %function lineCallback(l, src, event)
        %end
        
        %function shiftCoord(l, shift)
        %end
    end
end

function freeze(ax)
set(ax, {'XLimMode', 'YLimMode', 'ZLimMode'}, {'manual', 'manual', 'manual'});
end

function defreeze(ax)
set(ax, {'XLimMode', 'YLimMode', 'ZLimMode'}, {'auto', 'auto', 'auto'});
end


function addPointContextMenu(pnt, el)
    ax = pnt.Parent;
    m  = uicontextmenu('Parent', ax.Parent, 'UserData', pnt);
    pnt.UIContextMenu = m;
    uimenu('Parent', m, 'Label', 'Delete', 'Callback', @(src, event)removeData(src, event, el));
end

function ix = getPointInsertIx(l, p)
n = numel(l.XData);
if  n <= 1
    ix = n+1;
else
    v  = [l.XData - p(1), l.YData - p(2)];
    [v1, v2]   = deal(v(1:end-1,:), v(2:end, :));
    [nv1, nv2] = deal(sqrt(dot(v1,v1,2)), sqrt(dot(v2,v2,2)));
    costh = dot(v1, v2, 2)./(nv1.*nv2);
    [val, ix] = min(costh);
    ix = ix +1;
    if val > 0 % only accept if th > pi/2
        % choose first or last
        if nv1(1) < nv2(end)
            ix = 1;
        else
            ix = n+1;
        end
    end
end
end

function removeData(src, event, l)
set(src.Parent.UserData, {'XData', 'YData', 'ZData'}, {[], [], []});
l.cleanUpPoints();
l.resetLinePoints();
end



function flag = validData(l)
[nx, ny, nz] = deal(numel(l.XData), numel(l.YData), numel(l.ZData));
flag =  nx == ny && (isempty(nz) || ny == nz);
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
        
