classdef UIItem < handle
    properties
        Tag = '';
        panel
        fullHeight
        titleHeight
        layout
        fixedHeight = false;
        controls
        controlWidths
        dummyPanel = [];
        dummyStack = 'bottom';
        enableState
        collapseDirection = 'down';
        Children
    end
    properties (Dependent)
        Parent
        collapse
        fullTitle
        collapseTitle
        style
    end
    %% Appearance/Font props
    properties (Dependent)
        Title
        Position
        Visible
        Enable  = 'on';
        BackgroundColor
        ForegroundColor 
        titleColor
        titleBackgroundColor
        dummyColor
        FontName
        titleFontName
        FontSize
        titleFontSize
        FontWeight
        titleFontWeight
        FontAngle
        titleFontAngle
        HorizontalAlignment
    end
    
    %% Fixed props
    properties (SetAccess = immutable)
        Type      = 'UIItem';
        Units     = 'pixels';
        FontUnits = 'points'
    end
    %%
    methods
        %% constructor
        function d = UIItem(varargin)
            opt = struct('Parent',          [], ...
                         'controls',        {{{[]}}}, ...
                         'controlWidths',   {{}}, ...
                         'Callback',        '', ...
                         'Position', [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title', 'Empty UI-column item', ...
                         'style', 'default', ...
                         'reduceToUsefulHeight', true);
            
            [opt, extraOpt] = merge_options(opt, varargin{:}); % extra assumed to be props
            % currently requires at least one non-epmty child
            if all(cellfun(@isempty, opt.controls{1}))
                opt.controls{1}{1} = uicontrol('Style', 'text', 'String', char(9786*ones(1,10)));
            end
            
            % get parent
            if isempty(opt.Parent)
                p = gcf;
            else
                p = opt.Parent;
            end
            
            d.panel = uipanel('Parent', p, 'Units', 'pixels', 'Position', opt.Position, ...
                              'ButtonDownFcn', @d.toggleCollapse, 'Visible', 'off');
            
            % add cosmetic panel
            %if ~strcmp(d.style.name, 'simple')
            %    d.dummyPanel = uipanel('Parent', d.panel, 'Units', d.Units, 'Position', d.panel.Position - [0 0 0 d.titleHeight], ...
            %                           'BorderWidth', 0, 'HandleVisibility', 'off', 'Visible', d.panel.Visible);
            %end
            
            % assign controls and set layout
            d = assignControls(d, opt.controls);
            d.controlWidths = opt.controlWidths;
            
            d.fullHeight  = 0;
            
            d.setStyle(opt.style, extraOpt);
            d.Title       = opt.Title;
            d.Enable = 'all';
            
            d.panel.SizeChangedFcn = @d.updateLayout;
            d.Visible = 'on';
            
            % reduce height if requested, call after visible as inner
            % position tends to change
            if opt.reduceToUsefulHeight
                d.reduceHeight();
            end
            d.Visible = opt.Visible;
        end
    
        %% Appearance/Font props
        function set.Title(d, val)
            if ~d.collapse
                str = [char(9660), ' ', val];
            else
                str = [char(9658), ' ', val];
            end
            d.panel.Title = str;
        end
        function val = get.Title(d)
           val = d.panel.Title(3:end);
        end
        %------------------------------------------------------------------
        function set.Position(d,val)
            om = d.layout.params.outerMargins;
            newpos = val + [om(1), om(3), -(om(1)+om(2)), -(om(3)+om(4))];
            d.panel.Position = newpos;
            ipos = d.panel.InnerPosition;
            if ipos(end) == 0
                d.panel.SizeChangedFcn();
            end
            if ~isempty(d.dummyPanel)
                dom     = d.layout.params.dummyOuterMargins;
                if strcmp('top', d.dummyStack)
                    d.dummyPanel.Position = max(0, ...
                        [dom(1)+1, dom(3)+1, val(3)-(dom(1)+dom(2)), val(4)-(dom(3)+dom(4)+d.titleHeight)] );
                else
                    d.dummyPanel.Position = max(0, val + ...
                    [dom(1), dom(3), -(dom(1)+dom(2)), -(dom(3)+dom(4))] );
                end
            end
        end
%                 
%             else
%                 dm = d.layout.params.dummyMargins;
%                 if all(dm>=0) % dummy panel is in front (has panel as parent)
%                     d.dummyPanel.Position = ...
%                         [dm(1)+1, dm(3)+1, val(3)-(dm(1)+dm(2)),  max(0, val(4)-(dm(3)+dm(4)+d.titleHeight))];
%                     d.panel.Position = val;
%                 else          % % dummy panel is behind (has same parent as panel)
%                     [dmp, dmn] = deal(dm.*(dm>0), -(dm.*(dm<0)));
%                     d.dummyPanel.Position = val + ...
%                         [dmp(1), dmp(3), -(dmp(1)+dmp(2)),  -(dmp(3)+dmp(4))];
%                     d.panel.Position = val + ...
%                         [dmn(1), dmn(3), -(dmn(1)+dmn(2)),  -(dmn(3)+dmn(4))];
%                 end
%             end
%         end
        function val = get.Position(d)
            om  = d.layout.params.outerMargins;
            val = d.panel.Position - [om(1), om(3), -(om(1)+om(2)), -(om(3)+om(4))]; 
%             if isempty(d.dummyPanel) || all(d.layout.params.dummyMargins >= 0)
%                 val = d.panel.Position;
%             else
%                 [pp, dpp] = deal(d.panel.Position, d.dummyPanel.Position);
%                 val = [min(pp(1:2), dpp(1:2)), max(pp(3:4), dpp(3:4))];
%                %[pp; dpp]
%             end
        end
        %------------------------------------------------------------------
        function set.Visible(d, val)
            if any(strcmp(val, {'on', 'off'}))
                d.panel.Visible = val;
                if ~isempty(d.dummyPanel)
                    d.dummyPanel.Visible = val;
                end
                set(d.panel.Children, 'Visible', val);
            end
        end
        function val = get.Visible(d)
            vp = d.panel.Visible;
            if all(strcmp(vp, {d.panel.Children.Visible}))
                val = vp;
            else
                val = 'mixed';
            end
        end
        %------------------------------------------------------------------
        function set.Enable(d, val)
            if any(strcmp(val, {'on', 'off', 'all', 'inactive'}))
                % if previous was 'on', save to enableState
                 c = d.panel.Children;
%                 get(c(isprop(c, 'Enable')), 'Enable')
                if strcmp(val, 'off')
                    d.enableState = get(c, 'Enable');
                    %d.enableState = arrayfun(@(c)c.Enable, d.panel.Children, 'UniformOutput', false);
                end
                if strcmp(val, 'on') && ~isempty(d.enableState) % set to previous enableState
                    if ischar(d.enableState)
                        c.Enable = d.enableState;
                    else
                        for k = 1:numel(c)
                            c(k).Enable = d.enableState{k};
                        end
                    end
                else
                    if strcmp(val, 'all'), val = 'on'; end
                    d.setProp('Enable', val);
                   % set(d.panel.Children, 'Enable', val);
                end
            end
        end
        function val = get.Enable(d)
            val = d.getProp('Enable');
%             c = d.panel.Children;
%             if ~isempty(c)
%                 e1 = c(1).Enable;
%                 if all(strcmp(e1, {c.Enable}))
%                     val = e1;
%                 else
%                     val = 'mixed';
%                 end
%             end
        end
        %------------------------------------------------------------------
        function set.BackgroundColor(d, val)
            d.panel.BackgroundColor = val;
            %set(d.Children, 'BackgroundColor', val);
        end
        function val = get.BackgroundColor(d)
            val = d.getProp('BackgroundColor');
        end

        function set.ForegroundColor(d, val)
            %d.panel.ForegroundColor = val;
            set(d.Children, 'ForegroundColor', val);
        end
        function val = get.ForegroundColor(d)
            val = d.getProp('ForegroundColor');
        end
        
        function set.titleColor(d, val)
            d.panel.ForegroundColor = val;
        end
        function val = get.titleColor(d)
            val = d.panel.ForegroundColor;
        end
        
        function set.titleBackgroundColor(d, val)
            d.panel.BackgroundColor = val;
        end
        function val = get.titleBackgroundColor(d)
            val = d.panel.BackgroundColor;
        end
        
        function set.dummyColor(d, val)
            if ~isempty(d.dummyPanel)
                d.dummyPanel.BackgroundColor = val;
            end
        end
        function val = get.dummyColor(d)
            if ~isempty(d.dummyPanel)
                val = d.dummyPanel.BackgroundColor;
            else
                val = 'none';
            end
        end
        
        %------------------------------------------------------------------
        
        function set.FontName(d, val)
            set(d.Children, 'FontName', val);
        end
        function val = get.FontName(d)
            val = d.getProp('FontName');
        end
        
        function set.titleFontName(d, val)
            d.panel.FontName = val;
        end
        function val = get.titleFontName(d)
            val = d.panel.FontName;
        end
        
        function set.FontSize(d, val)
            set(d.Children, 'FontSize', val);
        end
        function val = get.FontSize(d)
            val = d.getProp('FontSize');
        end
        
        function set.titleFontSize(d, val)
            d.panel.FontSize = val;
        end
        function val = get.titleFontSize(d)
            val = d.panel.FontSize;
        end
        
        function set.FontWeight(d, val)
            set(d.Children, 'FontWeight', val);
        end
        function val = get.FontWeight(d)
            val = d.getProp('FontWeight');
        end
        
        function set.titleFontWeight(d, val)
            d.panel.FontWeight = val;
        end
        function val = get.titleFontWeight(d)
            val = d.panel.FontWeight;
        end
        
        function set.FontAngle(d, val)
            set(d.Children, 'FontAngle', val);
        end
        function val = get.FontAngle(d)
            val = d.getProp('FontAngle');
        end
        
        function set.titleFontAngle(d, val)
            d.panel.FontAngle = val;
        end
        function val = get.titleFontAngle(d)
            val = d.panel.FontAngle;
        end
        %------------------------------------------------------------------
        function set.HorizontalAlignment(d, val)
            d.setProp('HorizontalAlignment', val);
        end
        function val = get.HorizontalAlignment(d)
            val = d.getProp('HorizontalAlignment');
        end
        %------------------------------------------------------------------

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
                    
        %% Remaining set/get     
        function set.collapse(d, val)
            if ~val
                if strcmp(d.collapseDirection, 'down')
                    d.Position(4) = d.fullHeight;
                else
                    p = d.Position;
                    d.Position = [p(1), p(2)-(d.fullHeight-d.titleHeight), p(3), d.fullHeight];
                end
            else
                if strcmp(d.collapseDirection, 'down')
                    d.Position(4) = d.titleHeight;
                else
                    p = d.Position;
                    d.Position = [p(1), p(2)+(d.fullHeight-d.titleHeight), p(3), d.titleHeight];
                end
            end
            % also update title
            d.Title = d.Title;
        end
        function val = get.collapse(d)
            % may not be exacttly equal due to round-off 
            if abs(d.Position(4) - d.titleHeight) < .1
                val = true;
            else
                val = false;
            end
        end
        
        function set.Parent(d, val)
            d.panel.Parent = val;
        end
        function val = get.Parent(d)
            val = d.panel.Parent;
        end
        
        %% utility and callbacks (move some to static)
        function val = CurrentPoint(d, ~, ~) % relative to d.panel
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
               
        function toggleCollapse(d, ~, ~)
            p = d.CurrentPoint;
            if p(2) > d.Position(4) - d.titleHeight -1
                if d.collapse
                    d.collapse = false;
                else
                    d.collapse = true;
                end
            end
        end
        
        function updateLayout(d, ~, ~)
            if ~d.collapse
                % switch of SiceChangedFcn while positioning items 
                fn = d.panel.SizeChangedFcn;
                d.panel.SizeChangedFcn = '';
                d.fullHeight = d.Position(4);
                d = positionUIItems(d);
                d.panel.SizeChangedFcn = fn;
            end
            % also evaluate parents SizeChangedFcn if present
            %d.Position = d.panel.Position;
            if isprop(d.panel.Parent, 'SizeChangedFcn') && isa(d.panel.Parent.SizeChangedFcn, 'function_handle')
                feval(d.panel.Parent.SizeChangedFcn)
            end
        end
        
        function menu = listboxContextMenu(d, str, callbacks)
            % make default listbox context menu for uicontrol in figure f
            menu = uicontextmenu(figureParent(d));
            for k = 1:numel(str)
                uimenu(menu, 'Label', str{k}, 'Callback', callbacks{k});
            end
        end
        
%         function d = setPropertiesByLayout(d, opt)
%             switch d.style.name
%                 case 'simple'
%                     set(d, {'BackgroundColor', 'ForegroundColor', 'titleBackgroundColor', 'titleForegroundColor'}, ...
%                         {[0.94 0.94 0.94], [0.1328 0.2109 0.3496], [0.94 0.94 0.94],  [0.2656 0.4219 0.6992] });
%                     set(d, {'FontName', 'FontSize', 'titleFontSize', 'FontWeight', 'FontAngle'}, ...
%                         {'Helvetica', 10,        10,             'normal',     'normal'});
%                     d.HorizontalAlignment = 'left';
%                     
%                 case 'invert'
%                     set(d, {'BackgroundColor', 'ForegroundColor', 'titleBackgroundColor', 'titleForegroundColor'}, ...
%                         {[0.94 0.94 0.94], [0.1328 0.2109 0.3496], [0.2656 0.4219 0.6992], [0.94 0.94 0.94] });
%                     set(d, {'FontName', 'FontSize', 'titleFontSize', 'FontWeight', 'FontAngle'}, ...
%                         {'Helvetica', 10,        10,             'normal',     'normal'});
%                     d.HorizontalAlignment = 'left';     
%             end
%             for k = 1:2:numel(opt)
%                 d.(opt{k}) = opt{k+1};
%             end
%         end
        
        function d = reduceHeight(d) % remove extra y-space
            %nDummy = ~isempty(d.dummyPanel) && strcmp(d.dummyStack, 'top');
            pos = vertcat(d.Children.Position);
            om  = d.layout.params.outerMargins;
            m   = d.layout.params.margins;
            miny = min(pos(:,2));
            d.Position(4) = ceil(d.Position(4) - miny + sum(m(3:4)));%om(3) + om(4));
        end
        
        function d = setProp(d, prop, val)
            c = d.Children;
            set(c(isprop(c, prop)), prop, val);
        end
        function val = getProp(d, prop)
            c = d.Children;
            val = get(c(isprop(c, prop)), prop);
            if ~ischar(val) && numel(c) > 1
                if ~ischar(val{1}) || all(strcmp(val{1}, val(2:end)))
                    val = val{1};
                else
                    val = 'mixed';
                end
            end
        end
        
        function d = setTitleProp(d, prop, val)
            d.panel.(prop) = val;
        end
        function val = getTitleProp(d, prop)
            val = d.panel.(prop);
        end
        
        function d = assignControls(d, c) %set parent/children
            for k = 1:numel(c)
                for l = 1:numel(c{k})
                    if isprop(c{k}{l}, 'Parent')
                        c{k}{l}.Parent = d.panel;
                    end
                end
            end
            d.controls = c;
            % set Children (only panels of items)
            cld = horzcat(c{:});
            cld = cld(~cellfun(@isempty, cld));
            d.Children = vertcat(cld{:});
        end
        
        function set.style(d, styleName)
            d.setStyle(styleName);
        end
        function styleName = get.style(d)
            styleName = d.layout.styleName;
        end
        
        function d = setStyle(d, varargin)
            setUIStyle(d, varargin{:});
            %if ~isempty(d.panel.SizeChangedFcn)
            %    feval(d.panel.SizeChangedFcn);
            %end
            %d.Visible = 'on';
        end
        
        function d = positionUIItems(d)
            %disp(d.Title)
            %disp(d.Position)
            assert(isa(d.panel, 'matlab.ui.container.Panel'));
            
            params = d.layout.params;
            lineHeight      = ceil(params.lineHeightAt10*d.FontSize/10);
            %lineHeightPopup = 26*d.FontSize/10;
            lineHeightPopup = lineHeight +4;
            %d.titleHeight   = 19 + (d.titleFontSize - 10)*2;
            items = d.controls;
            nrows = numel(items);
            widths = d.controlWidths;
            if isempty(widths)
                widths = cell(1, nrows);
            end
            assert(numel(widths) == nrows);
            itemWidths     = cell(nrows, 1);
            rowHeights     = nan(nrows, 1);
            numLinesUseful = ones(nrows, 1);
            for k = 1:nrows
                l = items{k};
                rowItemWidths = nan(1, numel(l));
                curWidthAvailable = d.panel.Position(3) - (numel(l)-1)*params.hskip - sum(params.margins(1:2));
                input = widths{k};
                itemHeight = lineHeight*ones(numel(l), 1);
                if isempty(input), input = nan(1, numel(l)); end
                for m = 1:numel(l)
                    if isprop(l{m}, 'Type')
                        if strcmp(l{m}.Type, 'uicontrol')
                            switch l{m}.Style
                                case 'edit' % 'edit' uses fixed width
                                    rowItemWidths(m) = params.editWidth;
                                    %elseif strcmp(l{m}.Style, 'pushbutton')
                                    %    lineItemWidths(m) = layout.lineHeight;
                                case 'text'% use extent of text
                                    %rowItemWidths(m) = 60;
                                    l{m}.Extent(3);
                                case 'listbox'
                                    itemHeight(m) = nan;
                                    %numLinesUseful(m) = max(numLinesUseful(m), numel(l{m}.String));
                                case 'popupmenu' % popupmenu requires extra height
                                    itemHeight(m) = lineHeightPopup;
                                    %numLinesUseful(m) = max(numLinesUseful(m), numel(l{m}.String));
                                case 'pushbutton'
                                    rowItemWidths(m) = 1.6*lineHeight;
                            end
                        end
                    end
                    % reset width if given as input
                    if ~isnan(input(m))% inputWidths are relative
                        rowItemWidths(m) = curWidthAvailable*input(m);
                    end
                end
                % destribute remaining widths
                isn = isnan(rowItemWidths);
                if any(isn)
                    rowItemWidths(isn) = (curWidthAvailable - sum(rowItemWidths(~isn)))/nnz(isn);
                end
                itemWidths{k} = max(floor(rowItemWidths), 0);  % don't want negative
                % set height if fixed
                if ~any(isnan(itemHeight))
                    rowHeights(k) = max(itemHeight);
                end
            end
            % destribute remaining heights
            availableHeight  = d.panel.Position(4) - d.titleHeight ...
                                - (numel(items)-1)*params.vskip - sum(params.margins(3:4));
            isn = isnan(rowHeights);
            remainingHeight = availableHeight - sum(rowHeights(~isn));
            if any(isn)
                rowHeights(isn) = remainingHeight/nnz(isn);
            end
            rowHeights = max(floor(rowHeights), 0);
            
            % loop through again and position from top accoring to computed widths/heights
            sx = params.margins(1) + 1;
            curp = [sx, d.panel.Position(4) - d.titleHeight - params.margins(3) + params.vskip + 1];
            for k = 1:nrows
                curp = [sx, curp(2)-params.vskip-rowHeights(k)+1];
                for m = 1:numel(items{k})
                    if ~isempty(items{k}{m})
                        items{k}{m}.Position = ([curp, itemWidths{k}(m), rowHeights(k)]);
                        % some adjustments (set text closer to item below)
                        %if strcmp(items{k}{m}.Type, 'uicontrol')
                        %    if strcmp(items{k}{m}.Style, 'text')
                                %fzpix  = items{k}{m}.FontSize * 0.014 * 96; % approx fontsize in pixels
                                %adjust = max((lineHeight - fzpix)/2, 0);
                                %items{k}{m}.Position(2) = items{k}{m}.Position(2) - adjust;
                        %    end
                        %end
                    end
                    curp(1) = curp(1) + params.hskip + itemWidths{k}(m);
                end
            end
        end
    end
end



function param = layoutParameters()
param = struct(...
    'lineHeightAt10', 20, ...
    'vskip', 2, ...    % vertical space between items
    'hskip', 2, ...    % horizontal ...
    'lskip', 5, ...    % skip before leftmost item
    'rskip', 5, ...    % ... rightmost
    'tskip', 2, ...    % top
    'bskip', 2, ...    % bottom
    'panelWidth', 300, ... % just some default
    'tab', 0, ...
    'vPanelSkip', 5, ...
    'editWidth', 60);  % fixed
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
