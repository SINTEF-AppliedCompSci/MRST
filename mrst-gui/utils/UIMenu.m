classdef UIMenu < handle
    properties
        items 
        panel
        layout
        fullHeight
        width
        titleHeight = 18;
        collapseDirection
        %fixUpper          = true;
        level = 1;
        dummyPanel
        dummyStack = 'bottom';
        %lines
        Children
    end
    
    properties (Dependent)
        Parent
        collapse
        fullTitle
        collapseTitle
        style
    end
    
    properties (Dependent)
        Title
        Value = [];
        Position
        Visible
        Enable
        BackgroundColor
        ForegroundColor
        FontName
        FontSize
        FontWeight
        FontAngle
        %HorizontalAlignment
    end
    
    methods
        function m = UIMenu(varargin)
            opt = struct('Parent',          [], ...
                         'Callback',        'disp(''Hello'')', ...
                         'Position',  [10 10 200 400], ...
                         'Visible',         'on', ...
                         'Title', '', ...
                         'style', 'default', ...
                         'items', {{}}, ...
                         'itemFieldNames', {{}}, ...
                         'collapseDirection', 'up');
            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            if ~isempty(opt.Parent)
                p = opt.Parent;
            else
                p = gcf;
            end
                     
            %assert(strcmp(p.Type, 'figure'), 'Parent of uiMenu must be a figure');   
            
            m.panel  = uipanel('Parent',            p, ...
                               'Units',             'pixels', ...
                               'Position',          opt.Position, ...
                               'ButtonDownFcn',     @m.resize, ...
                               'SizeChangedFcn',    @m.destributeItems, ...
                               'Visible',           'off');
            
            m.items = opt.items;
            % assign items 
            for k = 1:numel(m.items)
                m.items{k}.panel.Parent     = m.panel;
                if isa(m.items{k}, 'UIMenu')
                    m.items{k}.level = m.level + 1;
                end
                % also point sub-panels to their item
                m.items{k}.panel.UserData = m.items{k};
            end
            % set children (main panels of items)
            cld = cellfun(@(x)x.panel, m.items, 'UniformOutput', false);
            m.Children = vertcat(cld{:});
            
            m.fullHeight  = 0; % trigger first draw
            
            %m.panel.Position(4);
            
            m.width = m.panel.Position(3);
            m.setStyle(opt.style, extraOpt);
            m.width = m.Position(3); 
            m.Title       = opt.Title;
            m.collapseDirection = opt.collapseDirection;
            if strcmp(m.dummyStack, 'top')
                m.dummyPanel.ButtonDownFcn = @m.resize;
            end
            m.Visible = 'on';
           
            
            
     
            %m.destributeItems();
        end
        
        %% Appearance/Font props
        function set.Title(m, val)
            if ~m.collapse
                str = [char(9660), ' ', val]; % char(9650), char(9651)
            else
                str = [char(9658), ' ', val]; % char(9660), char(9661)
            end
            %end
            m.panel.Title = str;
        end
        function val = get.Title(m)
            val = m.panel.Title(3:end);
        end
        %------------------------------------------------------------------
        function set.Position(d,val)
            if ~isempty(d.dummyPanel)
                dom     = d.layout.params.dummyOuterMargins;
                if strcmp('top', d.dummyStack)
                    d.dummyPanel.Position = max(-inf, ...
                        [dom(1), dom(3), val(3)-(dom(1)+dom(2)), max(0, val(4)-(dom(3)+dom(4)+d.titleHeight))] );
                else
                    d.dummyPanel.Position = max(-inf, val + ...
                        [dom(1), dom(3), -(dom(1)+dom(2)), -(dom(3)+dom(4))] );
                end
            end
            om = d.layout.params.outerMargins;
            d.panel.Position = val + [om(1), om(3), -(om(1)+om(2)), -(om(3)+om(4))];
            ipos = d.panel.InnerPosition;
            if ipos(end) == 0
                d.panel.SizeChangedFcn();
            end
            %             if isempty(d.dummyPanel)
            %                 d.panel.Position = val;
            %             else
            %                 dm = d.layout.params.dummyMargins;
            %                 if all(dm>=0) % dummy panel is in front (has panel as parent)
%                     d.dummyPanel.Position = ...
%                         [dm(1)+1, dm(3)+1, val(3)-(dm(1)+dm(2)),  max(0,val(4)-(dm(3)+dm(4)+d.titleHeight))];
%                     d.panel.Position = val;
%                 else          % % dummy panel is behind (has same parent as panel)
%                     [dmp, dmn] = deal(dm.*(dm>0), -(dm.*(dm<0)));
%                     d.dummyPanel.Position = val + ...
%                         [dmp(1), dmp(3), -(dmp(1)+dmp(2)),  -(dmp(3)+dmp(4))];
%                     d.panel.Position = val + ...
%                         [dmn(1), dmn(3), -(dmn(1)+dmn(2)),  -(dmn(3)+dmn(4))];
%                 end
%             end
        end
        function val = get.Position(d)
            om  = d.layout.params.outerMargins;
            val = d.panel.Position - [om(1), om(3), -(om(1)+om(2)), -(om(3)+om(4))];
            
            %             if isempty(d.dummyPanel) || all(d.layout.params.dummyMargins >= 0)
            %                 val = d.panel.Position;
            %             else
            %                 [pp, dpp] = deal(d.panel.Position, d.dummyPanel.Position);
            %                 val = [min(pp(1:2), dpp(1:2)), max(pp(3:4), dpp(3:4))];
            %             end
        end
        %------------------------------------------------------------------
        function set.Visible(m, val)
            if any(strcmp(val, {'on', 'off'}))
                m.panel.Visible = val;
                if ~isempty(m.dummyPanel)
                    m.dummyPanel.Visible = val;
                end
                setPropChildren(m.panel, 'Visible', val);
            end
        end
        function val = get.Visible(m)
            val = m.panel.Visible;
        end
        %------------------------------------------------------------------
        function set.Enable(m, val)
            if any(strcmp(val, {'on', 'off', 'all', 'inactive'}))
                % recursive, use set/get of items
                % if item is also of class uiMenu, Matlab does not allow call of
                % set-function, make call in seperate function (hack)
                for k = 1:numel(m.items)
                    setProp(m.items{k}, 'Enable', val);
                end
            end
        end
        function val = get.Enable(m)
            vals = cell(1, numel(m.items));
            for k = 1:numel(m.items)
                vals{k} = getProp(m.items{k}, 'Enable');
            end
            
            if all(strcmp(vals{1}, vals))
                val = vals{1};
            else
                val = 'mixed';
            end
        end
        %------------------------------------------------------------------        
        function set.Parent(m, val)
            m.panel.Parent = val;
        end
        function val = get.Parent(m)
            val = m.panel.Parent;
        end
        %------------------------------------------------------------------    
            
        function set.BackgroundColor(m, val)
            m.panel.BackgroundColor = val;
        end
        function val = get.BackgroundColor(m)
            val = m.panel.BackgroundColor;
        end
        
        function set.ForegroundColor(m, val)
            m.panel.ForegroundColor = val;
        end
        function val = get.ForegroundColor(m)
            val = m.panel.ForegroundColor;
        end
        
        function set.FontName(m, val)
            m.panel.FontName = val;
        end
        function val = get.FontName(m)
            val = m.panel.FontName;
        end
        
        function set.FontSize(m, val)
            m.panel.FontSize = val;
        end
        function val = get.FontSize(m)
            val = m.panel.FontSize;
        end
        
        function set.FontWeight(m, val)
            m.panel.FontWeight = val;
        end
        function val = get.FontWeight(m)
            val = m.panel.FontWeight;
        end
        
        function set.FontAngle(m, val)
            m.panel.FontAngle = val;
        end
        function val = get.FontAngle(m)
            val = m.panel.FontAngle;
        end
        
        function set(d, nm, val)
            if ischar(nm)
                nm  = {nm};
                val = {val};
            end
            
            for k = 1:numel(d)
                for l = 1:numel(nm)
                    d(k).(nm{l}) = val{l};
                end
            end
        end
        
        %------------------------------------------------------------------
        function set.collapse(m, val)
            if ~val
                if strcmp(m.collapseDirection, 'down')
                    m.Position(4) = m.fullHeight;
                else
                    p = m.Position;
                    m.Position = [p(1), p(2)-(m.fullHeight-m.titleHeight), p(3), m.fullHeight];
                end
            else
                if strcmp(m.collapseDirection, 'down')
                    m.Position(4) = m.titleHeight;
                else
                    p = m.Position;
                    m.Position = [p(1), p(2)+(m.fullHeight-m.titleHeight), p(3), m.titleHeight];
                end
            end
            % also update title
            m.Title = m.Title;
        end
        function val = get.collapse(m)
            % may not be exacttly equal due to round-off
            if (m.panel.Position(4) - m.titleHeight) < .1
                val = true;
            else
                val = false;
            end
        end
        
        function toggleCollapse(d, ~, ~)
            p = d.CurrentPoint;
            if p(2) > d.panel.Position(4) - d.titleHeight -1
                if d.collapse
                    d.collapse = false;
                else
                    d.collapse = true;
                end
            end
        end
        
        function setChildren(d, prop, val)
            fn = fieldnames(d.items);
            for k =1:numel(fn)
                d.items.(fn{k}).(prop) = val;
            end
        end
        
        function destributeItems(d, ~, ~)
            %only destribute items if x-extent has changed and if d is not collapsed
            if ~d.collapse
                %disp('.')
                %l = getGuiLayoutDefaults('layout');
                mrgs  = d.layout.params.margins;
                omrgs = d.layout.params.outerMargins;
                isItem = cellfun(@(x)isa(x, 'UIItem'), d.items);
                vskips = d.layout.params.vskipItem*isItem + ...
                    d.layout.params.vskipMenu*(~isItem);
                %vskip = d.layout.params.vskip;
                % only apply vskip below uiitems (not uimenus)
                %vskip_all = vskip * ones(numel(d.Children), 1);%arrayfun(@(x)isa(x.UserData, 'uiItem'), d.panel.Children);
                heights = cellfun(@(x)x.Position(4), d.items);
                d.fullHeight = d.titleHeight + sum(heights) + sum(vskips) + sum(mrgs(3:4)) + 0*sum(omrgs(3:4));
                %d.fullHeight = d.getFullHeight(vskip_all);
                %if ~d.collapse
                dh = d.Position(4) - d.fullHeight;
                %else
                %    dh = d.Position(4) - d.titleHeight;
                %end
                if (dh ~= 0) || (d.width ~= d.Position(3)) || ...
                    (d.level > 1 && (d.Position(3) ~= d.Parent.Position(3) - sum(mrgs(1:2))))
                    % switch of SizeChangedFcn while drawing to prevent calling from
                    % items
                    %disp(d.Title)
                    %disp(d.Position)
                    
                    fn = d.panel.SizeChangedFcn;
                    d.panel.SizeChangedFcn = '';
                    
                    d.Position(4) = d.Position(4) - dh;
                    d.width      = d.Position(3);
                    if d.level == 1
                        d.Position(2) = d.Position(2) + dh;
                        %else
                        % d.Position(3) = d.width;
                    end
                    
                    %c = d.Children;
                    % last half of children may be dumm panels - dont loop over
                    % them
                    psize     = d.Position(3:4);
                    itemWidth = psize(1)-sum(omrgs(1:2))-sum(mrgs(1:2));
                    %dm = d.layout.params.dummyMargins(1:2);
                    %dm = dm.*(dm<0);
                    cury  = psize(2)-0*omrgs(4)-mrgs(4)-d.titleHeight+0*vskips(1)+1;% + dm(2);
                    for k = 1:numel(d.items) %numel(d.items):-1:1
                        % use item Postition, not item.panel
                        %item = c(k).UserData;
                        %csize = d.items{k}.Position(3:4);
                        if k > 1
                            vskip = vskips(k-1);
                        else
                            vskip = 0;
                        end
                        cury = cury - heights(k) - vskip;
                        d.items{k}.Position = [mrgs(1)+1, cury, itemWidth, heights(k)];
                    end
                    d.panel.SizeChangedFcn = fn;
                    %                [c(1).Position; c(2).Position];
                end
            end
            % update parent
            if ~(d.level == 1) && isprop(d.panel.Parent, 'SizeChangedFcn') && isa(d.panel.Parent.SizeChangedFcn, 'function_handle')
                feval(d.panel.Parent.SizeChangedFcn)
            end
            
        end
        
        %         function h = getFullHeight(d, vskip_all)
        %             if numel(d.items) > 0
        %                 heights = cellfun(@(x)x.Position(4), d.items);
        %                 %dm = d.layout.params.dummyMargins;
        %                 h = d.titleHeight + sum(heights) + sum(vskip_all) + sum(d.layout.params.margins(3:4));
        %             else
        %                 h = d.titleHeight;
        %             end
        %         end
        
        function val = CurrentPoint(d, ~, ~)
            %val = d.panel.Parent.CurrentPoint - d.panel.Position(1:2);
            p = d.panel;
            d = [0 0];
            while ~isa(p, 'matlab.ui.Figure')
                d = d + p.Position(1:2);
                p = p.Parent;
            end
            val = p.CurrentPoint - d;
        end
        
        function f = figureParent(d)
            f = d.panel.Parent;
            while ~isa(f, 'matlab.ui.Figure')
                f = f.Parent;
            end
        end
        
        function set.style(m, styleName)
            m.setStyle(styleName);
        end
        function styleName = get.style(m)
            styleName = m.layout.styleName;
        end
        
        function m = setStyle(m, varargin)
            setUIStyle(m, varargin{:})
            for k = 1:numel(m.items)
                m.items{k}.setStyle(varargin{:});
            end
            m.destributeItems()
        end
        
        function m = insertItem(m, item, pos)
            % insert menu item at specified position (or at end if defaulted)
            ni = numel(m.items);
            if nargin < 3
                pos = ni+1;
            end
            assert(pos <= ni+1, 'Menu item position exceeds total number of items')
            item.Parent = m.panel;
            item.level  = m.level +1;
            m.items     = [m.items(1:pos-1), {item}, m.items(pos:end)];
            m.Children  = [m.Children(1:pos-1); item.panel; m.Children(pos:end)];
            m.destributeItems();
        end
             
        function resize(m, ~, ~)
            p = m.CurrentPoint;
            f = figureParent(m);
            motionFcn     = f.WindowButtonMotionFcn;
            upFcn         = f.WindowButtonUpFcn;
%             if false %p(1) >= d.panel.Position(3)-5
%                 offset = d.panel.Position(3)-p(1);
%                 f.WindowButtonMotionFcn = {@dragHorisontally, offset};
%                 f.WindowButtonUpFcn     = {@clearMotionFcn, motionFcn};
            if m.Position(4) - p(2) < m.titleHeight
                m.toggleCollapse();
            else % check if between items
                dy = arrayfun(@(x)x.Position(2), m.Children)-p(2);
                ix = and(dy > 0,dy <= 10);
                if any(ix)  % resize item above
                    %check that current item is not collapsed
                    item = m.Children(ix).UserData;
                    if ~isa(item, 'UIItem')
                        return;
                    end
                    if ~item.collapse
                        set(f,'Pointer','top')
                        offset = m.Children(ix).Position(2)-p(2);
                        f.WindowButtonMotionFcn = @dragVertically;
                        f.WindowButtonUpFcn     = @clearMotionFcn;
                    end
                end
            end
                    
            function dragHorisontally(src, event, offset)
                p = m.CurrentPoint+offset;
                m.panel.Position(3) = max(p(1), 20);
            end
            
            function dragVertically(src, event)
                p = m.CurrentPoint+offset;
                newHeight = item.Position(4) + ...
                           (item.Position(2) - p(2));
                item.Position(4) = max(newHeight, item.titleHeight+10);
                %item.panel.Position(2) = max(p(2), 10);
                %feval(d.panel.SizeChangedFcn);
            end
            
            function clearMotionFcn(src, event)
                f.WindowButtonMotionFcn = motionFcn;
                f.WindowButtonUpFcn     = upFcn;
                set(f,'Pointer','arrow')
            end
        end
        
        %function clickThroughDummy(m, ~, ~)
            
    end
end



function setPropChildren(obj, prop, val)
if isprop(obj, 'Children')
    c = obj.Children;
    for k = 1:numel(c)
        if isprop(c(k), prop)
            c(k).(prop) = val;
            setPropChildren(c(k), prop, val);
        end
    end
end
end
  
function setProp(obj, prop, val)
obj.(prop) = val;
end

function val = getProp(obj, prop)
val = obj.(prop);
end

function s = uniquify(s)
    [su, ia, ib] = unique(s);
    if numel(su) ~= numel(s)
        for k = 1:numel(su)
            ix = find(ia(k)==ib);
            if numel(ix) > 1
                for l = 1:numel(ix)
                    s{ix(l)} = [s{ix(l)}, num2str(l)];
                end
            end
        end
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
