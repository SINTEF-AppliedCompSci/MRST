classdef InteractiveLineOnSurf < handle
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
        cutDir    = [0 1 0];
        colorData = [];
        subset    = [];
        slice     = [];
        G       
        sliceType = 'along';  % 'normalPlane'
    end
     properties (Dependent)
        XData
        YData
        ZData
     end
    
    methods
        function l = InteractiveLineOnSurf(G, varargin)
            l.G = G;
            opt       = struct('Parent',   [], ...
                               'Callback', @(src, event)'Hello', ...
                               'colorData', []);
            [opt, other] = merge_options(opt, varargin{:});
            
            pnts = struct('XData', [], 'YData', [], 'ZData', []);
            [pnts, other] = merge_options(pnts, other{:});
            l.Callback  = opt.Callback;
            l.colorData = opt.colorData;
            
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
                                 'Parent', l.Parent, 'PickableParts', 'none', ...
                                 'Color', 'm', 'LineWidth', 4);
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
                markers = [];
            end
            l.markers = markers;
            
            l.resetLinePoints();
%            l.updateCSys();
            l.addAxisContextMenu();
            
            l.trajectory.Visible = 'on';
            set(l.markers, 'Visible', 'on');
            l.updateSlice();
        end
        
        function val = get.XData(l)
            if ~isempty(l.markers)
                val = vertcat(l.markers.XData);
            else
                val = [];
            end
        end
        function set.XData(l, val)
            if ~(numel(val) == numel(l.markers))
                error('Got %d values for %d markers', numel(val), numel(l.markers));
            else
                for k = 1:numel(l.markers)
                    l.markers(k).XData = val(k);
                end
            end
        end
                    
        
        function val = get.YData(l)
            if ~isempty(l.markers)
                val = vertcat(l.markers.YData);
            else
                val = [];
            end
        end
        function set.YData(l, val)
            if ~(numel(val) == numel(l.markers))
                error('Got %d values for %d markers', numel(val), numel(l.markers));
            else
                for k = 1:numel(l.markers)
                    l.markers(k).YData = val(k);
                end
            end
        end
        
        function val = get.ZData(l)
            if ~isempty(l.markers)
                val = vertcat(l.markers.ZData);
            else
                val = [];
            end
        end
        function set.ZData(l, val)
            if ~(numel(val) == numel(l.markers))
                error('Got %d values for %d markers', numel(val), numel(l.markers));
            else
                for k = 1:numel(l.markers)
                    l.markers(k).ZData = val(k);
                end
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
                                     'ButtonDownFcn', @l.movePoint, 'Parent', l.Parent);
            addPointContextMenu(h, l)
            h.MarkerFaceColor = 'r';
            if ~strcmp(l.drawMode, 'along')
                l.markers(min(n, nd)).MarkerFaceColor = 'g';
            end
            l.markers = [l.markers(1:n-1), h, l.markers(n:nd)];
            l.resetLinePoints();
            l.updateSlice();
        end
        
        function deletePoint(l, n)
            nd =  numel(l.XData);
            if nargin < 2
                n =nd;
            end
            assert(n>=1 && n <= nd);
            
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
                cp = l.cursorPointOnSurf();
                if ~any(isnan(cp))
                    %cp = ax.CurrentPoint(1,1:2);
                    lp = [src.XData, src.YData, src.ZData];
                    motionFcn = fg.WindowButtonMotionFcn;
                    upFcn     = fg.WindowButtonUpFcn;
                    offset = lp-cp;
                    fg.WindowButtonMotionFcn = {@dragPoint, src, offset};
                    fg.WindowButtonUpFcn     = {@clearMotionFcn, motionFcn, upFcn};
                end
            end
            function dragPoint(src, event, pnt, offset)
                %ax.CurrentPoint
                %newp = ax.CurrentPoint(1, 1:2)+offset;
                newp = l.cursorPointOnSurf(offset);
                if ~any(isnan(newp))
                    [xl, yl, zl] = deal(ax.XLim, ax.YLim, ax.ZLim);
                    set(pnt, {'XData', 'YData', 'ZData'}, ...
                        {max(xl(1), min(newp(1), xl(2))), max(yl(1), min(newp(2), yl(2))), ...
                        max(zl(1), min(newp(3), zl(2)))});
                    l.resetLinePoints;
                end
            end
            
             function clearMotionFcn(src, ~, motionFcn, upFcn)
                src.WindowButtonMotionFcn = motionFcn;
                src.WindowButtonUpFcn     = upFcn; 
                %l.resetLinePoints;
                l.updateSlice();
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
                l.Callback([], []);
            else
                %l.delete();
            end
        end
        
        function redistributePoints(l, nPnts)
            npOld = numel(l.XData);
            if npOld >= 2 && npOld ~= nPnts
                t  = l.trajectory;
                npt = numel(t.XData);
                p = interp1((0:npt-1)'/(npt-1), [t.XData(:), t.YData(:), t.ZData(:)], (0:nPnts-1)/(nPnts-1));
                for k = 1:max(nPnts, npOld)
                    if k <= nPnts
                        if k <= npOld
                            set(l.markers(k), {'XData', 'YData', 'ZData'}, ...
                                {p(k,1), p(k,2), p(k,3)});
                        else
                            l.insertPoint(p(k,:));
                            l.markers(k).MarkerFaceColor = 'g';   
                        end
                    else
                        l.deletePoint();
                    end
                end
                l.resetLinePoints();
            end
        end
               
        function addAxisContextMenu(l)
            ax = l.Parent;
            if isprop(ax, 'ContextMenu')
                cmenu = 'ContextMenu';
            elseif isprop(ax(1), 'UIContextMenu')
                cmenu = 'UIContextMenu';
            end
            m = ax.(cmenu);
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
        end
        
        function putPoint(l, src, event)
            %ax = l.Parent;
            %p = ax.CurrentPoint(1,1:2);
            %p = event.IntersectionPoint;
            p = l.cursorPointOnSurf();
            if event.Button == 1 && ~any(isnan(p))
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
            if ~isempty(l.slice) && isvalid(l.slice)
                delete(l.slice);
            end
            if isvalid(l.Parent)
                ax = l.Parent;
                if isprop(ax, 'ContextMenu')
                    cmenu = 'ContextMenu';
                elseif isprop(ax, 'UIContextMenu')
                    cmenu = 'UIContextMenu';
                end
                if numel(ax.(cmenu).Children) >= 2
                    delete(ax.(cmenu).Children(1));
                else
                    l.Parent.(cmenu) = gobjects(0);
                end
                l.Parent.ButtonDownFcn = '';
                
                if isvalid(l.Parent.Parent)
                    l.Parent.Parent.WindowButtonMotionFcn = '';
                    l.Parent.Parent.WindowButtonUpFcn = '';
                    l.Parent.Parent.Pointer = 'arrow';
                end
            end
        end
        
        function p = cursorPointOnSurf(l, offset)
            if nargin < 2
                offset = [0 0 0];
            end
            % find intersection of current point view-line and surface
            % defined by trajectory and cutDir. Equation is
            %  p = s1 + x1*ds = ti + x2*dti + x3*dc
            % s is view-segment, t is trajectory, dc is cutDir
            p = nan(1,3);
            seg  = l.Parent.CurrentPoint + ones(2,1)*offset;
            traj = [l.trajectory.XData(:), l.trajectory.YData(:), ...
                l.trajectory.ZData(:)];
            % find curve segment for intersection by projecting to null-space
            % of view-segment and cutDir
            dp = diff(traj);
            ds = diff(seg);
            dc = l.cutDir;
            if abs(ds*dc') < .99*norm(ds)*norm(dc)
                n  = null([dc;ds]);
                x2  = ( bsxfun(@minus, seg(1,:), traj(1:end-1,:))*n )./(dp*n);
                % need 0 <= t <= 1
                ix = x2 >= 0 & x2 <= 1;
                % extrapolate
                if ~isempty(x2)
                    ix(1)   = x2(1) <= 1;
                    ix(end) = x2(end) >= 0;
                end
                ix = find(ix);
                if numel(ix) >= 2 % select closest to viewer
                    pi = (1-x2(ix))'*traj(ix,:) + x2(ix)'*traj(ix+1,:);
                    dv = bsxfun(@minus, pi, seg(1,:));
                    [~, ii] = min(sum(dv.^2,2));
                    ix = ix(ii);
                end
                if ~isempty(ix)
                    x2 = x2(ix);
                    % project view-dir onto null-space of cutDir
                    n = (ds - (dc*ds')*dc)';
                    x1 = (([1-x2, x2]*traj(ix:ix+1,:)-seg(1,:))*n)/(diff(seg)*n);
                    p  = [1-x1, x1]*seg;
                end
            end
        end
                
        function updateSlice(l)
            c = getTrajPoints(l);
            poly = computeGridSlicePolygons(l.G, c, 'cutDir', l.cutDir);
            if isempty(poly)
                %warning('Trajectory does not intersect grid');
                return
            end
            if isempty(l.colorData)
                if isempty(l.slice)
                    l.slice = patch('Faces', poly.nodes, 'Vertices', poly.coords3D, ...
                        'FaceColor', 'y', 'FaceAlpha', .2, 'PickableParts', 'none', ...
                        'Parent', l.Parent);
                else
                    set(l.slice, {'Faces', 'Vertices'}, {poly.nodes, poly.coords3D});
                end
            else
                if isempty(l.slice)
                    l.slice = patch('Faces', poly.nodes, 'Vertices', poly.coords3D, ...
                        'FaceColor', 'flat', 'FaceAlpha', .2, 'PickableParts', 'none', ...
                        'Parent', l.Parent, 'FaceVertexCData', l.colorData(poly.cellIx));
                else
                    set(l.slice, {'Faces', 'Vertices', 'FaceVertexCData'}, ...
                        {poly.nodes, poly.coords3D, l.colorData(poly.cellIx)});
                end
            end
        end
        
        function setCutDir(l, d)
            l.cutDir = d;
            l.updateSlice();
        end
        
        function setAngle(l, v)
            p  = [l.XData([1,end]), l.YData([1,end]), l.ZData([1,end])];
            n  = diff(p);
            n  = n/norm(n);
            v1 = [1 0 0];
            if abs(n(3)) > 10*max(abs(n(1:2))) % approx vertical
                v2 = [0 1 0];
            else
                if abs(n(2)) < abs(n(1))
                    v1 = [0 1 0];
                end
                v1(1:2) = v1(1:2) - (n(1:2)*v1(1:2)')*n(1:2);
                v1 = v1/norm(v1);
                v2 = null([v1;n])';
            end
            l.setCutDir(cos(deg2rad(v))*v1 + sin(deg2rad(v))*v2);
        end
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



function c = getTrajPoints(l, extrap)
if nargin < 2
    extrap = 1;
end
atol = 1;
c = [l.trajectory.XData(:), l.trajectory.YData(:), l.trajectory.ZData(:)];
% project c onto plane given by p = mean(c), n = cutDir
[p, n] = deal(mean(c), l.cutDir);
pp = c - (bsxfun(@minus, c, p)*n')*n;
% remove parts pointing in opposite direction
md = diff(pp([1,end],:));
ok = true(size(pp,1),1);
kcur = 1;
for k = 2:numel(ok)
    if dot(md, diff(pp([kcur,k],:))) > 0
        kcur = k;
    else
        ok(k) = false;
    end
end
pp = pp(ok,:);
% resample
len = [0; cumsum(sqrt(sum(diff(pp).^2, 2)))];
c  = samplePiecewiseLinear(@(x)interp1(len/len(end), pp, x), [0 1], atol);

if extrap > 0
    lims = get(l.Parent, {'XLim', 'YLim', 'ZLim'});
    lims = vertcat(lims{:});
    db   = extrap*sqrt(sum(diff(lims, [], 2).^2));
    % start
    ds = diff(c(1:2, :));
    c(1,:) = c(1,:) - ds*(db/norm(ds));
    % end 
    de = diff(c(end-1:end,:));
    c(end,:) = c(end,:) + de*(db/norm(ds));
end

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
        
