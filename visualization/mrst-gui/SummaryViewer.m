classdef SummaryViewer < handle
    % SYNOPSIS:
    %   d = SummaryViewer()                     : select one or more summary files
    %   d = SummaryViewer({fn1, fn2, ...})      : input name(s) of summary files
    %   d = SummaryViewer({smry1, smry2, ...})  : input summary structures as obtained
    %                                             by readEclipseSummaryUnFmt
    %
    % DESCRIPTION:
    %   Reads (optionally) and displays data from one or more summary file(s) in GUI  
    % 
    % Customization:
    % Some properties are not exposed in the GUI, but can be set by 
    % command. Assuming the the viewer has been started with
    % >> d = SummaryViewer(...)
    % these include:
    % d.LineWidth   : line width for all plots (deafult 2)
    % d.FontSize    : font size used in axis (deafult 10)
    % d.MarkerSize  : size of any marker used in plot (default 8)
    % d.xDataType   : x-data used in plot. Valid options are 
    %                 'date' (default) and 'step' (i.e., ministep number) 
    %
    % Order of colors, line-styles and markers can be set respectively in 
    % properties colorOrder, lineStyleOrder and markerOrder. For instance
    % >> d.colorOrder = hsv(20);
    % resets line-colors and takes effect in the next interaction with the
    % GUI or immediatly by typing
    % >> d.update();
    
    properties
        Figure
        SubItemFigure
        Axes
        selector
        optionSelector
        smry
        casePaths
        layout = struct('menuWidth', 240, 'subItemWidth', 150)
        markerOrder = {'none', 'x', 'o', 's', 'd'};
        lineStyleOrder = {'-', '--', '-.', ':'};
        colorOrder
        subItemIx
        subItemSelectors
        currentSubItems = {};
        keepAxesLims  = true;
        filterZero    = false;
        lineOptions   = struct('LineWidth', 2, ...
                               'MarkerSize', 8, ...
                               'xDataType',  'date');
    end
    properties (Dependent)
        caseNames 
        caseSelection
        orderCasesBy 
        orderNamesBy 
        orderPropsBy 
        selectSubItems
        scaleValues
        LineWidth
        MarkerSize
        FontSize
        xDataType
    end
    
    methods
        function d = SummaryViewer(varargin)
            [d.smry, d.caseNames, d.casePaths]= handleInput(varargin);
            if isempty(d.smry{1}), return; end
            d.caseSelection = (1:numel(d.smry))';
            d.Figure = figure('Name', 'main');
            d.Axes   = axes('Parent', d.Figure, 'Units', 'pixels');
            if isempty(d.colorOrder)
                d.colorOrder = d.Axes.ColorOrder;
            end
            d.selector = SummarySelectorExtended(d.smry, 'Parent', d.Figure, 'Units', 'pixels');
            d.optionSelector = SummaryOptionsSelector('caseNames', d.caseNames, 'Parent', d.Figure, 'Units', 'pixels');
            d.Figure.SizeChangedFcn = @d.updateLayout;
            d.Figure.CloseRequestFcn = @d.closeFigure;
            d.selector.Callback = @d.update;
            d.optionSelector.Callback = @d.update; 
            d.updateLayout();
            box(d.Axes, 'on')
            grid(d.Axes, 'on')
            d.Axes = addAxesContextMenu(d.Axes);
        end
        
        %------------------------------------------------------------------
        function set.caseNames(d, nameList)
            d.optionSelector.caseNames = nameList;
        end
        function nameList = get.caseNames(d)
            nameList = d.optionSelector.caseNames;
        end  
        function set.caseSelection(d, val)
            d.optionSelector.caseSelection = val;
        end
        function val = get.caseSelection(d)
            val = d.optionSelector.caseSelection;
        end 
        %------------------------------------------------------------------
        function set.orderCasesBy(d, str)
            d.optionSelector.lineProperties{1} = str;
        end
        function str = get.orderCasesBy(d)
             str = d.optionSelector.lineProperties{1};
        end
        function set.orderNamesBy(d, str)
            d.optionSelector.lineProperties{2} = str;
        end
        function str = get.orderNamesBy(d)
             str = d.optionSelector.lineProperties{2};
        end
        function set.orderPropsBy(d, str)
            d.optionSelector.lineProperties{3} = str;
        end
        function str = get.orderPropsBy(d)
             str = d.optionSelector.lineProperties{3};
        end
        %------------------------------------------------------------------
        function set.selectSubItems(d, flag)
            d.optionSelector.doSelectSubItems = flag;
        end
        function flag = get.selectSubItems(d)
            flag = d.optionSelector.doSelectSubItems;
        end
        function set.scaleValues(d, flag)
            d.optionSelector.doScale = flag;
        end
        function flag = get.scaleValues(d)
            flag = d.optionSelector.doScale;
        end
        %------------------------------------------------------------------
        function set.LineWidth(d, val)
            d.lineOptions.LineWidth = val;
            ix = isprop(d.Axes.Children, 'LineWidth');
            set(d.Axes.Children(ix), 'LineWidth', val);
        end
        function val = get.LineWidth(d)
            val = d.lineOptions.LineWidth;
        end
        function set.MarkerSize(d, val)
            d.lineOptions.MarkerSize = val;
            ix = isprop(d.Axes.Children, 'MarkerSize');
            set(d.Axes.Children(ix), 'MarkerSize', val);
        end
        function val = get.MarkerSize(d)
            val = d.lineOptions.MarkerSize;
        end
        function set.FontSize(d, val)
            d.Axes.FontSize = val;
        end
        function val = get.FontSize(d)
            val = d.Axes.FontSize;
        end
        function set.xDataType(d, str)
            if any(strcmpi(str, {'date', 'step'}))
                d.lineOptions.xDataType = str;
                d.update();
            end
        end
        function str = get.xDataType(d)
            str = d.lineOptions.xDataType;
        end
        
        %------------------------------------------------------------------
        function update(d, src, event)
            if nargin < 2
                src.Tag = '';
                event = [];
            end
            [nms, prps] = deal(d.selector.curNames, d.selector.curProps);
            ax = d.Axes;
            if d.keepAxesLims
                cla(ax);
            else
                cla(ax, 'reset');
                box(ax, 'on')
                grid(ax, 'on')
            end
            if d.selectSubItems && ~isempty(d.SubItemFigure) && isvalid(d.SubItemFigure) && ~strcmp(src.Tag, 'subItemSelector')
                d.clearSubItemSelectors();
            elseif ~d.selectSubItems && ~isempty(d.SubItemFigure) && isvalid(d.SubItemFigure)
                close(d.SubItemFigure);
            end
            count   = 0;
            caseCnt = 1;
            [data, info, propIx] = deal(cell(nnz(d.caseSelection)*numel(nms)*numel(prps), 1));
            for kc = 1:numel(d.smry)
                if d.caseSelection(kc)
                    nameCnt = 1;
                    caseInc = 0;
                    for kn = 1:numel(nms)
                        propCnt = 1;
                        nameInc = 0;
                        for kp = 1:numel(prps)
                            if any(strcmp(prps{kp}, d.smry{kc}.getKws(nms{kn})))
                                [caseInc, nameInc] = deal(1);
                                [vals, unit, spec] =  d.smry{kc}.get(nms{kn}, prps{kp}, ':', true);
                                isMultiple = size(vals, 1) > 1;
                                if isempty(vals)
                                    fprintf('No valid data for %s : %s\n', nms{kn},  prps{kp});
                                    continue; 
                                end
                                count = count + 1;
                                if d.filterZero
                                    vals(vals==0) = nan;
                                end
                                if ~d.selectSubItems || size(vals,1) == 1
                                    ix = 1:size(vals,1);
                                    propIx{count} = repmat([propCnt, nameCnt, caseCnt], [numel(ix), 1]);
                                    propCnt = propCnt +1;
                                else
                                    if ~isempty(spec)
                                        d.updateSubItems(src, event, nms{kn}, prps{kp}, spec);
                                        ix = lookupHelper(d.subItemIx, nms{kn}, prps{kp}, 'vals');
                                        if max(ix) > size(vals,1)
                                            fprintf('Non-compatible index selection (likely output mismatch between cases)\n');
                                            ix = ix(ix <= size(vals,1));
                                        end
                                    else
                                        ix = 1;
                                    end
                                    propIx{count} = [((1:numel(ix))' + propCnt -1) ,repmat([nameCnt, caseCnt], [numel(ix), 1]) ];
                                    propCnt = propCnt + numel(ix);
                                end
                                vals = vals(ix,:);
                                sc = 1;
                                if d.scaleValues
                                    sc = max(max(abs(vals)));
                                    if sc == 0, sc = 1;end
                                end
                                data{count} = vals/sc;
                                info{count} = getPlotInfo(d, kc, nms{kn}, prps{kp}, unit, ix, sc, isMultiple);
                            end
                        end
                        nameCnt = nameCnt + nameInc;
                    end
                    caseCnt = caseCnt + caseInc;
                end
            end
            % plotting
            if count > 0
                info = info(1:count);
                caseNo = cellfun(@(s)s.caseNo, info);
                h = cell(numel(d.smry), 1);
                for kc = 1:numel(d.smry)
                    ii = caseNo == kc;
                    if any(ii)
                        if strcmpi(d.lineOptions.xDataType, 'date')
                            xdata = d.selector.time{kc};
                        else
                            xdata = d.smry{kc}.ministep;
                        end
                        h{kc} = line(ax, xdata, vertcat(data{ii}), 'LineWidth', d.LineWidth, 'MarkerSize', d.MarkerSize);
                    end
                end
                h = vertcat(h{:});
                propIx = vertcat(propIx{1:count});
                ix = getLinePropIx(d, numel(prps), numel(prps), numel(d.smry), propIx);
                [mo, lso] = deal(d.markerOrder(:), d.lineStyleOrder(:));
                set(h, {'Marker'},    mo(ix{1}), ...
                       {'Color'},     num2cell(d.colorOrder(ix{2},:), 2), ...
                       {'LineStyle'}, lso(ix{3}) );
                [lstr, ltitle, ix] = getLegendInput(d, info, 40);
                legend(h(ix), lstr, 'Interpreter', 'none', 'FontName', 'FixedWidth', 'FontWeight', 'bold');
                title(d.Axes.Legend, ltitle);
                xtickangle(ax, 30)
                if numel(d.selector.curProps) == 1 && strcmp(d.selector.curProps{1}, 'WMCTL')
                    d.Axes.YLim = [-1, 8];
                    d.Axes.YTick = (0:7);
                    ax.YTickLabels = {'SHUT/STOP', 'ORAT', 'WRAT', 'GRAT', 'LRAT', 'RESV', 'THP', 'BHP'};
                    ytickangle(ax, 60)
                    d.keepAxesLims = false;
                else
                    d.Axes.YTickMode = 'auto';
                    d.Axes.YTickLabelMode = 'auto';
                    d.Axes.YLimMode = 'auto';
                    ytickangle(ax, 0)
                    d.keepAxesLims = true;
                end
            end
        end
        
        %------------------------------------------------------------------
        function updateLayout(d, varargin)
            fpos = d.Figure.Position;
            mw = d.layout.menuWidth;
            if d.optionSelector.collapse
                optH = d.optionSelector.titleHeight;
            else
                optH = numel(d.caseNames)*(d.optionSelector.titleHeight+6) + ...
                       5.5*d.optionSelector.titleHeight;
            end
            d.optionSelector.Position = [5, 10, mw, optH];
            d.selector.Position = max(0, [5, optH + 15, mw, fpos(4)-optH - 20]);
            d.Axes.Position     = max(0, [mw + 50, 50, fpos(3)-mw-75, fpos(4)-75]);
        end

        %------------------------------------------------------------------        
       
        
%         function updateOptions(d, varargin)
%             d.caseNames      = d.optionSelector.caseNames;
%             d.caseSelection  = d.optionSelector.caseSelection;
%             d.selectSubItems = d.optionSelector.doSelectSubItems;
%             d.scaleValues    = d.optionSelector.doScale;
%             [d.orderCasesBy, d.orderNamesBy, d.orderPropsBy] = ...
%                 deal(d.optionSelector.lineProperties{:});
%             if nargin > 1
%                 d.update(varargin{:});
%             end
%         end
        
        %------------------------------------------------------------------        
        function updateSubItems(d, src, event, nm, prp, info)
            if isempty(d.SubItemFigure) || ~isvalid(d.SubItemFigure)
                d.SubItemFigure = limitedToolbarFigure;
                fpos = d.Figure.Position;
                d.SubItemFigure.Position(1:2) = [fpos(1) + fpos(3)/2, fpos(2)]; 
                d.SubItemFigure.CloseRequestFcn = @d.closeSubItemsFigure;
                d.SubItemFigure.SizeChangedFcn = @d.updateSubItemLayout;
            end
            [nm_, prp_] = deal(matlab.lang.makeValidName(nm), matlab.lang.makeValidName(prp));
            if ~isfield(d.subItemIx, nm_) || ~isfield(d.subItemIx.(nm_), prp_)
                d.subItemIx.(nm_).(prp_) = struct('active', false, 'vals', 1, 'name', nm, 'prop', prp);
            end
            if ~d.subItemIx.(nm_).(prp_).active
                d.addSubItemSelector(src, event, nm, prp, nm_, prp_, info);
                d.subItemIx.(nm_).(prp_).active = true;
            end
        end
        
        %------------------------------------------------------------------        
        function addSubItemSelector(d, src, event,  nm, prp, nm_, prp_, info)
            str = getItemListStr(prp, info);
            nsel = numel(d.subItemSelectors);
            cur  = nsel +1;
            d.subItemSelectors{cur} = TimeStepSelector('Parent', d.SubItemFigure, 'tsteps', str, ...
                'Title', [nm, '-', prp]);
            d.subItemSelectors{cur}.Value = d.subItemIx.(nm_).(prp_).vals;
            d.updateSubItemLayout(src, event, true);
            d.currentSubItems = [d.currentSubItems; {nm_, prp_}];
            d.subItemSelectors{cur}.Callback = @(src, event)d.updateSubItemSelection(src, event, cur);
            d.subItemSelectors{cur}.panel.ButtonDownFcn = '';
            % reset context menu text
            d.subItemSelectors{cur}.selector.ContextMenu.Children(1).Text = 'Clear all';
            d.subItemSelectors{cur}.selector.ContextMenu.Children(2).Text = 'Select all';
            d.subItemSelectors{cur}.selector.Tag = 'subItemSelector';
        end
        
        %------------------------------------------------------------------        
        function updateSubItemLayout(d, src, event, setFigWidth)
            if nargin < 4
                setFigWidth  =false;
            end
            fpos = d.SubItemFigure.Position;
            mw = d.layout.subItemWidth;
            if setFigWidth
                ncol = max(2, min(numel(d.subItemSelectors), 10.5));
                d.SubItemFigure.Position(3) = (mw+2)*ncol;
            end
            mw = d.layout.subItemWidth;
            for k = 1:numel(d.subItemSelectors)
                d.subItemSelectors{k}.Position = [(k-1)*(mw+2), 1, mw, fpos(4)];
            end
        end
        
        %------------------------------------------------------------------
        function updateSubItemSelection(d, src, event, num)
            d.subItemIx.(d.currentSubItems{num,1}).(d.currentSubItems{num,2}).vals = d.subItemSelectors{num}.Value;
            d.update(src, event);
        end
        
        %------------------------------------------------------------------
        function clearSubItemSelectors(d, src, event)
            for k = 1:size(d.currentSubItems,1)
                d.subItemIx.(d.currentSubItems{k,1}).(d.currentSubItems{k,2}).active = false;
            end
            clf(d.SubItemFigure);
            [d.currentSubItems, d.subItemSelectors] = deal({});
        end
        
        %------------------------------------------------------------------
        function closeFigure(d, src, event)
            if ~isempty(d.SubItemFigure)
                delete(d.SubItemFigure);
            end
            delete(d.Figure);
        end
        
        %------------------------------------------------------------------
        function closeSubItemsFigure(d, src, event)
            d.clearSubItemSelectors();
            delete(d.SubItemFigure);
            d.optionSelector.doSelectSubItems = false;
            d.update(src, event);
        end
        
        %------------------------------------------------------------------
        function toggleSubItems(d, src, event)
            set(src.Parent.Children, 'Checked', false)
            src.Checked = true;
            d.selectSubItems = strcmp(src.Text, 'On');
            d.update(src, event);
        end
        
        function toggleScale(d, src, event)
            set(src.Parent.Children, 'Checked', false)
            src.Checked = true;
            d.scaleValues = strcmp(src.Text, 'On');
            d.update(src, event);   
        end
        
        function switchLineProps(d, src, event, fld, prp)
            set(src.Parent.Children, 'Checked', false)
            src.Checked = true;
            switch fld
                case 'case'
                    d.orderCasesBy = prp;
                case 'name'
                    d.orderNamesBy = prp;
                case 'prop'
                    d.orderPropsBy = prp;                    
            end
            d.update(src, event);
        end
        
        function setLineWidth(d, src, event)
            d.lineWidth = str2double(src.Text);
            d.update(src, event);
        end
        
        function switchFontSize(d, src, event)
            d.Axes.FontSize = str2double(src.Text);
        end
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [smry, names, paths] = handleInput(varargin)
if nargin == 0 || isempty(varargin{1})
    smry = selectSummaryFiles();
else
    smry = varargin{1}{1};
end
if ~iscell(smry)
    smry = {smry};
end
[names, paths] = deal([]);
hasFolderNames = false;
if ischar(smry{1}) && strcmp(smry{1}, 'folder')
    [smry, names] = selectFolder();
    hasFolderNames = true;
end
if all(cellfun(@ischar, smry))
    paths = smry;
end
for k = 1:numel(smry)
    if ischar(smry{k})
        [pth, nm] = fileparts(smry{k});
        if ~hasFolderNames
            names{k} = nm; %#ok
        end
        smry{k} = readEclipseSummaryUnFmt(fullfile(pth, nm));
    end
end
if isempty(names)
    names = applyFunction(@(k)sprintf('Case %k', k), (1:numel(smry)));
end
opt = struct('caseNames', {names});
opt = merge_options(opt, varargin{2:end});
names = opt.caseNames;
end

%--------------------------------------------------------------------------
function inp = selectSummaryFiles()  
inp = [];
[fnl, pthl] = uigetfile('*.SMSPEC', 'Select ECLIPSE summary specification(s)  (SMSPEC)', 'Multiselect', 'on');
    if ~isnumeric(fnl)
        if ~iscell(fnl), fnl = {fnl}; end
        inp = cell(1, numel(fnl));
        pth = pthl;
        for k = 1:numel(fnl)
            if iscell(pthl)
                pth = pthl{k};
            end
            inp{k} = fullfile(pth, fnl{k});
        end
    end
end

%--------------------------------------------------------------------------
function v = lookupHelper(s, f1, f2, fld)
v = s.(matlab.lang.makeValidName(f1)).(matlab.lang.makeValidName(f2)).(fld);
end

%--------------------------------------------------------------------------
function str = getItemListStr(prp, info)
[n, ni] = size(info);
ndig = ceil(log10(max(info)))+1;
switch prp(1)
    case {'B', 'C'} % block-format
        ss = @(ijk)sprintf('[%*d, %*d, %*d]', ndig(1), ijk(1), ndig(2), ijk(2), ndig(3), ijk(3));
    case 'R'
        if size(info,2) == 1
            ss = @(i)sprintf('[%*d]', ndig(1), i);
        elseif size(info,2) == 2
            ss = @(rr)sprintf('[%*d -> %*d]', ndig(1), rr(1), ndig(2), rr(2));
        end 
    otherwise
        if ni >= 1
            ss = @(i)sprintf('[%*d]', ndig(1), i);
        else
            ss = @(i)sprintf('[%*d]', 4, i);
        end
end
str = cell(1, n);
nIx = ceil(log10(n))+3;
for k = 1:n
    str{k} = sprintf('%*d : %s', nIx, k, ss(info(k,:)));
end
end

%--------------------------------------------------------------------------
function ix = getLinePropIx(d, nprp, nnms, ncs, k)
lProps = {'Marker', 'Color', 'LineStyle'};
props  = {d.orderPropsBy, d.orderNamesBy, d.orderCasesBy};
n      = [nprp, nnms, ncs];
ix = cell(1,3);
nlp = [numel(d.markerOrder), size(d.colorOrder,1), numel(d.lineStyleOrder)];
for j = 1:3
    isPrp = matches(props, lProps{j});
    prd   = [1, cumprod((n(1:2)-1).*isPrp(1:2) + 1)];
    fac   = isPrp.*prd;
    ix{j} = mod((k-1)*fac(:), nlp(j)) +1;
end
end

%--------------------------------------------------------------------------
function info = getPlotInfo(d, caseNo, name, prop, unit, specNo, scale, isMultiple)
info.caseNo = caseNo;
info.numLines = numel(specNo);
info.case =  d.caseNames(caseNo);
info.name = {name};
info.prop = {prop};
if scale ~= 1
    unit = sprintf('%7.2E %s', scale, unit);
end
info.unit = {sprintf('[%s]', unit)};    
if d.selectSubItems && isMultiple
    n = numel(specNo);
    info.case = repmat(info.case, [n,1]);
    info.name = repmat(info.name, [n,1]);
    info.unit = repmat(info.unit, [n,1]);
    info.prop = applyFunction(@(k)sprintf('%s(%*d)', prop, ceil(log10(max(specNo)+.5)), k), specNo(:));
end
end

%--------------------------------------------------------------------------
function [lstr, title, ix] = getLegendInput(d, info, maxNum)
nl = cellfun(@(s)s.numLines, info);
if d.selectSubItems
    ix = 1:sum(nl);
else
    ix = [1; 1+cumsum(nl(1:end-1))];
end
cn = applyFunction(@(s)s.case, info);
nm = applyFunction(@(s)s.name, info);
pr = applyFunction(@(s)s.prop, info);
un = applyFunction(@(s)s.unit, info);
M = {vertcat(cn{:}), vertcat(nm{:}), vertcat(pr{:}), vertcat(un{:})};
for k = 1:numel(M)
    nchar = max(cellfun(@numel, M{k}));
    M{k}  = applyFunction(@(s)sprintf('%*s', nchar, s), M{k});
end
M{3} = applyFunction(@(s1,s2)sprintf('%s %s', s1, s2), M{3}, M{4});
isSingle = cellfun(@(x)numel(unique(x))==1, M(1:3));
M = horzcat(M{1:3});
if all(isSingle), isSingle = ~isSingle; end
if all(isSingle) || ~any(isSingle)
    title = '';
    pat = '%s %s %s';
elseif nnz(isSingle) == 2
    title= sprintf('%s %s', M{1, isSingle});
    pat = '%s';
else
    title = sprintf('%s', M{1, isSingle});
    pat = '%s %s';
end
lstr = applyFunction(@(k)sprintf(pat, M{k, ~isSingle}), (1:size(M,1)));
if maxNum < numel(ix)
    fprintf('Displaying first %d legend entries\n', maxNum)
    lstr = lstr(1:maxNum);
    ix   = ix(1:maxNum); 
end
end

% -----------------------------------------------------------------
function ax = addAxesContextMenu(ax)
m = uicontextmenu('Parent', ax.Parent);
uimenu(m, 'Label', 'Export to new figure', 'Callback', @(src, event)copyToFig(src, event, ax))
ax.UIContextMenu = m;
end

% -----------------------------------------------------------------
function copyToFig(src, event, ax) %#ok
    f_new = figure;
    new = copyobj([ax, ax.Legend], f_new);
    ax  = findobj(new, 'Type', 'axes');
    drawnow
    set(ax, 'Units', 'normalized');
    set(ax, 'OuterPosition', [0 0 1 1],'Position',[.13 .11 .775 .815])
    set(get(gca,'Children'),'HitTest','on');
end

% -----------------------------------------------------------------
function [inp, names] = selectFolder()
[inp, names] = deal([]);
pth = uigetdir(pwd,'Select parent folder for (multiple) summary cases');
if ~isnumeric(pth)
    sub = dir(fullfile(pth, '*', '*.SMSPEC'));
    if ~isempty(sub)
        names = applyFunction(@(x)fullfile(filepartsname(x.folder), filepartsname(x.name)), sub);
        inp = applyFunction(@(x)fullfile(x.folder, x.name), sub);
        [names, inp] = deal(reshape(names, 1, []), reshape(inp, 1, []));
    end
    if numel(sub)>5
        ix = listdlg('ListString', names, 'ListSize', [300, 800]);%[300, max(250, 20*numel(sub))]);
        [names, inp] = deal(names(ix), inp(ix));
    end
end
end

% -----------------------------------------------------------------
function nm = filepartsname(pth)
    [~, nm] = fileparts(pth);
end

% -----------------------------------------------------------------



%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
