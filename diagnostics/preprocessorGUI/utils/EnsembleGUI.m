classdef EnsembleGUI < handle
    properties
        Figure
        Axes
        Menu
        Data
        m
        currentDiagnostics
        names
        injectors
        validMembers
        memberSelection     % zero for non-selected
        memNoToLocalNo      % mapping from ensemble number to subset number
        markers             % markers for cross plot
        labels              % marker labels
        lines               % lines for line-plots
        hists               % histograms for static variables
        fitLine             % best fitting line for crossplots
        crossPlotSelections = {'Model No', ...
            'Lorenz coefficient', ...
            'Sweep', ...
            'Recovery Factor', ...
            'NPV per volume', ...
            'maxNPV'};
        linePlotSelections  = {'F-Phi', ...
            'Sweep vs time', ...
            'Residence time distribution', ...
            'Oil Rates', ...
            'Water Rates', ...
            'WCUT vs time', ...
            'Recovery vs time', ...
            'Recovery factor vs Time', ...
            'NPV vs time'};
        histPlotSelections = {'Porosity','Perm.x', 'Perm.y', 'Perm.z'};
        layoutParams = struct('menuWidth',  300, ...
            'itemSpace', 40, 'includeToolTips', true);
        plotOpts       = [];
        injectorNames
        producerNames
        currentPlotType = 'cross';
        textColor
        highlightColors = get(groot, 'DefaultAxesColorOrder');
        diagnosticsViewerHandle = [];
    end
    properties (Dependent)
        xData
        yData
        xSelectionCrossPlot
        ySelectionCrossPlot
        selectionLinePlot
        selectionHistPlot
        selectionHist2Plot
        injectorIx
        producerIx
        sizeProp
        maxTOF
        timeUnit
        fluidProps
        npvProps
        useOriginalFluid
    end
    
    methods
        function d = EnsembleGUI(ensemble, varargin)
            opt = struct('memberIx', [], 'diagnosticsNo', 1, 'style', 'default');
            opt = merge_options(opt, varargin{:});
            ix = opt.diagnosticsNo;
            assert(isa(ensemble, 'ModelEnsemble'));
            if isempty(ensemble.diagnostics)
                warning('Ensemble has no computed diagnostics. Diagnostics must be precomputed by running ensemble.computeDiagnostics, returning ...');
                return
            end
            
            d.m = ensemble;
            % check for valid data
            d.validMembers = getValidMembers(d.m, ix, opt.memberIx);
            d.memberSelection = zeros(numel(d.validMembers), 1);
            d.memNoToLocalNo = zeros(d.m.nMembers,1);
            d.memNoToLocalNo(d.validMembers) = 1;
            d.memNoToLocalNo = cumsum(d.memNoToLocalNo);
            if isempty(d.validMembers)
                warning('Empty diagnostics for required realizations, returning ...');
                return
            end
            
            % Figure
            d.Figure = limitedToolbarFigure();
            
            % setup axes
            d.Axes = axes('Parent', d.Figure, 'Units', 'pixels');
            d.Axes = addAxesContextMenu(d.Axes);
            % name members according to
            d.names = arrayfun(@(i)[d.m.name, sprintf('%3.3d', i)], d.validMembers, ...
                'UniformOutput', false);
            
            % setup selectors/menu
            itemOpts = {'Parent', d.Figure, 'Visible','off', 'style', opt.style};
            msel = InputModelSelector('modelNames', d.names, 'Title', 'Highlight models', ...
                'ButtonString', 'Diagnostics Viewer', itemOpts{:});
            
            WP_sample = d.m.diagnostics(ix).WP{d.validMembers(1)};
            d.injectorNames    = {WP_sample.inj.name};
            d.producerNames   = {WP_sample.prod.name};
            wsel = WellSelector('injectors', d.injectorNames, 'producers', d.producerNames, ...
                'includeThreshold', false, itemOpts{:});
            
            csel = CrossplotSelector('props', {d.crossPlotSelections}, itemOpts{:}, ...
                'includeFreezeSwitch', true, 'includeSizePopup', true, ...
                'sizeProps', {{'none', 'allocation', 'volume'}});
            
            lsel = SimplePopupSelector('props', {d.linePlotSelections}, itemOpts{:});
            
            statistics = d.m.diagnostics(ix).statistics{1};
            Properties = struct('static',  struct('name', {d.histPlotSelections}, 'limits', ...
                {{[statistics.poro.globalmin,    statistics.poro.globalmax],...
                convertTo([statistics.perm.globalmin(1),  statistics.perm.globalmax(1)], milli*darcy),...
                convertTo([statistics.perm.globalmin(2),  statistics.perm.globalmax(2)], milli*darcy),...
                convertTo([statistics.perm.globalmin(3),  statistics.perm.globalmax(3)], milli*darcy)}}));
            
            
            hsel = PropertyHistogramSelector( itemOpts{:},...
                'props', Properties,...
                'includeFilter', true,...
                'Title','Histograms',...
                'includeLogSwitch', true);            
            
            tsel = ThresholdTOFSelector(itemOpts{:});
            fsel = TwoPhaseFluidSelector(itemOpts{:});
            nsel = TwoPhaseNPVSelector(itemOpts{:});
            
            m1   = UIMenu('Title', 'Settings', 'Parent', [], 'items', {fsel, nsel});
            
            d.Menu = UIMenu('Title', 'Menu', 'Parent', d.Figure, itemOpts{:}, ...
                'items', {msel, wsel, csel, lsel, tsel, m1, hsel});
            
            tt = [251 180 76; 58 152 216; 42 187 155; ...
                  252 121 122; 76 150 180; 102 128 50; 191, 79, 81]./255;
            tt2 = tatarizeMap;
            tt = [tt; tt2(3:end,:)];
            for k =1:numel(d.Menu.items)
                d.Menu.items{k}.BackgroundColor = tt(k,:);
                try
                    d.Menu.items{k}.titleColor = 'w';
                catch
                    d.Menu.items{k}.ForegroundColor = 'w';
                end
            end
            
            d.Figure.SizeChangedFcn        = @d.layout;
            d.layout();
            
            % process data
            d.Data = loadAndProcessData(d.m, ix, d.validMembers);
            
            wsel.communicationMatrix = d.Data.wellCommunication;
            wsel.Enable = 'on';
            
            d.emptyCurrentDiagnostics();
            d.crossPlotCallback();
            
            % setup callback functions
            msel.Callback = @d.memberSelectCallback;
            msel.modelViewerButtonCallback = @d.modelViewerButtonCallback;
            wsel.Callback = @d.regionSelectCallback;
            csel.Callback = @d.crossPlotCallback;
            lsel.Callback = @d.linePlotCallback;            
            tsel.Callback = @d.thresholdCallback;
            fsel.Callback = @d.fluidCallback;
            nsel.Callback = @d.npvCallback;
            hsel.Callback = @d.histPlotCallback;
            
            % setup markers
            d = setupPlotOpts(d);
            
            d.markers  = arrayfun(@(x)line('Visible', 'on', 'XData', [], 'YData', []), 1:numel(d.validMembers));
            d.labels   = arrayfun(@(x)text('Visible', 'off', 'Position', [.5 .5], ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold'), ...
                1:numel(d.validMembers));
            d.fitLine  = line('Visible', 'off', 'XData', [], 'YData', [], 'Color', [.6 .6 1], 'LineWidth', 2, 'LineStyle', '--');
            d.lines    = arrayfun(@(x)line('Visible', 'on', 'XData', [], 'YData', []), 1:numel(d.validMembers));
            d.hists    = arrayfun(@(x)line('Visible', 'on', 'XData', [], 'YData', []), 1:numel(d.validMembers));
            
            
            set(d.markers, d.plotOpts.point.prop, d.plotOpts.point.val);
            set(d.lines, d.plotOpts.line.prop, d.plotOpts.line.val);
            set(d.hists, d.plotOpts.line.prop, d.plotOpts.line.val);            % TO DO Properties for d.hists should be define later.
            
            
            set([d.markers, d.labels, d.lines, d.hists], 'ButtonDownFcn', @d.memberSelectCallback);
            for k = 1:numel(d.validMembers)
                rnum = d.validMembers(k);
                d.markers(k).Tag      = d.names{k};
                d.markers(k).UserData = rnum;
                d.labels(k).Tag       = d.names{k};
                d.labels(k).UserData  = rnum;
                d.labels(k).String    = num2str(rnum);
                d.lines(k).Tag        = d.names{k};
                d.lines(k).UserData   = rnum;
                d.hists(k).Tag        = d.names{k};
                d.hists(k).UserData   = rnum;
            end
            
            % setup text-colors
            d.textColor = struct('active',    csel.ForegroundColor, ...
                'inactive', [.75 .75 .75]);
            
            % initial selection
            wsel.autoCheckBox.Value = 1;
            wsel.injectorIx = (1:numel(d.injectorNames));
            wsel.producerIx = (1:numel(d.producerNames));
            lsel.ForegroundColor = d.textColor.inactive;
            
            % get useful limits of tof-thresholding
            rtd  = d.getDiagnosticsValues('RTD');
            pvi2year = mean(cellfun(@(x)x.volumes/(x.allocations*year), rtd));
            % set year-limits based on minimal max rtd-time
            tt = max(cellfun(@(r)max(r.t), rtd))/year;
            % but don't go above 50 pvi
            tt = min(tt, 10*pvi2year);
            
            tsel.yearLimits = [.1*pvi2year, tt];
            tsel.pviLimits  = [.1, tt/pvi2year];
            
            csel.xIx = 2;
            csel.yIx = 3;
            tsel.Value = 1;
            fsel.useOriginal = true;
            if d.useOriginalFluid
                d.computePhaseRates();
            end
            for i=2:numel(d.Menu.items)
                d.Menu.items{i}.collapse=true;
            end
        end
        
        %% get - functions
        function val = get.xData(d)
            val = cell2mat(get(d.markers, 'XData'));
        end
        
        function val = get.yData(d)
            val = cell2mat(get(d.markers, 'YData'));
        end
        
        function val = get.xSelectionCrossPlot(d)
            val = d.Menu.items{3}.xSelection;
        end
        
        function val = get.ySelectionCrossPlot(d)
            val = d.Menu.items{3}.ySelection;
        end
        
        function val = get.selectionLinePlot(d)
            val = d.Menu.items{4}.selection;
        end
        
        function val = get.selectionHistPlot(d)
            val.selection = d.Menu.items{7}.selection;
            val.logscale  = d.Menu.items{7}.logSwitchBox;
        end
        
        function val = get.injectorIx(d)
            val = d.Menu.items{2}.injectorIx;
        end
        
        function val = get.producerIx(d)
            val = d.Menu.items{2}.producerIx;
        end
        
        function val = get.sizeProp(d)
            val = d.Menu.items{3}.szSelection;
        end
        
        function val = get.maxTOF(d)
            val = d.Menu.items{5}.Value;
        end
        
        function val = get.timeUnit(d)
            val = d.Menu.items{5}.timeUnit;
        end
        
        function val = get.useOriginalFluid(d)
            val = d.Menu.items{6}.items{1}.useOriginal;
        end
        
        function val = get.fluidProps(d)
            fs = d.Menu.items{6}.items{1};
            val = struct('muw', fs.muw, 'muo', fs.muo, ...
                'swc', fs.swc, 'soc', fs.soc, ...
                'n',   fs.n);
        end
        
        function val = get.npvProps(d)
            ns = d.Menu.items{6}.items{2};
            val = struct('ro', ns.ro, 'rwp', ns.rwp, ...
                'rwi', ns.rwi, 'd', ns.d);
        end
        
        %% callback functions ---------------------------------------------
        function memberSelectCallback(d, src, event)                       %#ok<*INUSD>
            % handle select/deselect of ensemble members
            [addIx, removeIx] = deal([]);
            if nargin > 1 && any(strcmp(src.Type, {'line', 'text'}))
                memNo = d.memNoToLocalNo(src.UserData);
                if d.memberSelection(memNo) > 0 % remove selection
                    removeIx = memNo;
                else
                    addIx = memNo;
                end
            else % menu selection
                curmem = find(d.memberSelection);
                removeIx = setdiff(curmem, d.Menu.items{1}.ix);
                addIx    = setdiff(d.Menu.items{1}.ix, curmem);
            end
            if ~isempty(removeIx)
                d.memberSelection(removeIx) = 0;
            end
            if ~isempty(addIx)
                addVals = max(d.memberSelection) + (1:numel(addIx))';
                d.memberSelection(addIx) = addVals;
                uistack(d.markers(addIx), 'top');
                uistack(d.lines(addIx), 'top');
                uistack(d.labels(addIx), 'top');
            end
            [selIx, ~, colIx] = find(d.memberSelection);
            d.Menu.items{1}.ix = selIx;
            nCol = size(d.highlightColors, 1);
            if strcmp(d.currentPlotType, 'cross')
                [sp, sv] = deal(d.plotOpts.pointSelect.prop, d.plotOpts.pointSelect.val);
                [dp, dv] = deal(d.plotOpts.pointDeselect.prop, d.plotOpts.pointDeselect.val);
                
                set(d.markers(d.memberSelection==0), dp, dv);
                set(d.labels(d.memberSelection==0), 'Visible', 'off');
                set(d.markers(d.memberSelection>0), sp, sv);
                for k =1:numel(selIx)
                    col = d.highlightColors(mod(colIx(k), nCol)+1,:);
                    d.markers(selIx(k)).MarkerFaceColor = col;
                    if sum(col) > 1.5
                        d.labels(selIx(k)).Color = [.05 .05 .05];
                    else
                        d.labels(selIx(k)).Color = [.95 .95 .95];
                    end
                end
                set(d.labels(d.memberSelection>0), 'Visible', 'on');
            else
                [sp, sv] = deal(d.plotOpts.lineSelect.prop, d.plotOpts.lineSelect.val);
                [dp, dv] = deal(d.plotOpts.lineDeselect.prop, d.plotOpts.lineDeselect.val);
                set(d.lines(d.memberSelection>0), sp, sv);
                set(d.hists(d.memberSelection>0), sp, sv);
                for k =1:numel(selIx)
                    d.lines(selIx(k)).Color = d.highlightColors(mod(colIx(k), nCol)+1,:);
                    d.hists(selIx(k)).Color = d.highlightColors(mod(colIx(k), nCol)+1,:);
                end
                set(d.lines(d.memberSelection==0), dp, dv);
                set(d.hists(d.memberSelection==0), dp, dv);
            end
        end
        
        %------------------------------------------------------------------
        
        function crossPlotCallback(d, src, event)
            % handle cross-plot
            if ~strcmp(d.currentPlotType, 'cross')
                d.currentPlotType = 'cross';
                d.Menu.items{3}.ForegroundColor = d.textColor.active;
                d.Menu.items{4}.ForegroundColor = d.textColor.inactive;
                d.Menu.items{5}.ForegroundColor = d.textColor.inactive;
                d.Menu.items{6}.ForegroundColor = d.textColor.active;
                d.Menu.items{6}.items{1}.ForegroundColor  = d.textColor.active;
                d.Menu.items{6}.items{2}.ForegroundColor  = d.textColor.active;
                d.Menu.items{7}.ForegroundColor = d.textColor.inactive;
                set(d.lines, 'Visible', 'off')
                set(d.hists, 'Visible', 'off')
                d.Axes.XScale = 'linear';
                d.memberSelectCallback();
            end
            if nargin == 1 || strcmp(src.Tag, 'resize')
                if strcmp(d.sizeProp, 'none')
                    set(d.markers, 'MarkerSize', 14);
                else
                    % resize markers according to selected value
                    val = d.getDiagnosticsValues(d.sizeProp);
                    sz  = min(20, max(10, 15*round(val/mean(val))));
                    sz(d.memberSelection>0) = sz(d.memberSelection>0) + 4;
                    for k = 1:numel(d.markers)
                        d.markers(k).MarkerSize = sz(k);
                    end
                end
            end
            d.fitLine.Visible = 'off';
            % keep y-axes fixed when x-axis shows modelNo for better
            % viz/interactivity for showing development over time
            if strcmp(d.xSelectionCrossPlot, 'Model No')
                d.Axes.XLim = [0, d.m.nMembers+1];
                if ~any(strcmp(d.ySelectionCrossPlot, {'Model No', 'NPV per volume', 'maxNPV'}))
                    d.Axes.YLim = [0 1.05];
                elseif any(strcmp(d.ySelectionCrossPlot, {'Model No', 'NPV per volume'}))
                    d.Axes.YLimMode = 'auto';
                else
                    mx = d.getDiagnosticsValues('maxNPV');
                    mx = max([mx{:}, 1]);
                    d.Axes.YLim = [-mx mx];
                end
            else
                d.Axes.XLimMode = 'auto';
                d.Axes.YLimMode = 'auto';
            end
            
            set(d.markers, 'Visible', 'off')
            xSel = d.xSelectionCrossPlot;
            
            % handle freeze-options
            if ~d.Menu.items{3}.xFreeze
                x = d.getDiagnosticsValues(xSel);
                for k = 1:numel(d.markers)
                    d.markers(k).XData = x{k};
                    if ~isempty(x{k})
                        d.labels(k).Position(1) = x{k};
                    end
                end
                d.Axes.XLabel.String = xSel;
            elseif nargin > 1 && strcmp(src.Tag,'xFreeze')
                if numel(d.producerIx) == numel(d.producerNames)
                    pno = 'all';
                else
                    pno = regexprep(num2str(d.producerIx), '\s*', ', ');
                end
                if numel(d.injectorIx) == numel(d.injectorNames)
                    ino = 'all';
                else
                    ino = regexprep(num2str(d.injectorIx), '\s*', ', ');
                end
                d.Axes.XLabel.String = [xSel ,', injectors ', ino, ', producers ', pno];
            end
            ySel = d.ySelectionCrossPlot;
            if ~d.Menu.items{3}.yFreeze
                y = d.getDiagnosticsValues(ySel);
                for k = 1:numel(d.markers)
                    d.markers(k).YData = y{k};
                    if ~isempty(y{k})
                        d.labels(k).Position(2) = y{k};
                    end
                end
                d.Axes.YLabel.String = ySel;
            elseif nargin > 1 && strcmp(src.Tag,'yFreeze')
                if numel(d.producerIx) == numel(d.producerNames)
                    pno = 'all';
                else
                    pno = regexprep(num2str(d.producerIx), '\s*', ', ');
                end
                if numel(d.injectorIx) == numel(d.injectorNames)
                    ino = 'all';
                else
                    ino = regexprep(num2str(d.injectorIx), '\s*', ', ');
                end
                d.Axes.YLabel.String = [ySel ,', injectors ', ino, ', producers ', pno];
            end
            set(d.markers, 'Visible', 'on')
            
            % finally fit line if options is selected
            if d.Menu.items{3}.fitLine
                [x,y] = deal(d.xData, d.yData);
                ix = ~isnan(x) & ~isnan(y);
                if nnz(ix) > 1
                    [d.fitLine.XData, d.fitLine.YData] = ...
                        fitLineSegment(x(ix), y(ix), d.Axes.XLim, d.Axes.YLim);
                    d.fitLine.Visible = 'on';
                end
            end
        end
        
        % -----------------------------------------------------------------
        
        function linePlotCallback(d, src, event)
            % handle line plot
            d.fitLine.Visible = 'off';
            if ~strcmp(d.currentPlotType, 'line')
                d.currentPlotType = 'line';
                d.Menu.items{3}.ForegroundColor = d.textColor.inactive;
                d.Menu.items{4}.ForegroundColor = d.textColor.active;
                d.Menu.items{5}.ForegroundColor = d.textColor.active;
                d.Menu.items{6}.ForegroundColor = d.textColor.active;
                d.Menu.items{6}.items{1}.ForegroundColor  = d.textColor.active;
                d.Menu.items{6}.items{2}.ForegroundColor  = d.textColor.active;
                d.Menu.items{7}.ForegroundColor = d.textColor.inactive;
                set(d.markers, 'Visible', 'off')
                set(d.labels, 'Visible', 'off')
                set(d.hists, 'Visible', 'off')
                d.Axes.XScale = 'linear';
                d.memberSelectCallback();
            end
            set(d.lines, 'Visible', 'off')
            sel = d.selectionLinePlot;
            val = d.getDiagnosticsValues(sel);
            if numel(val) == 3
                if strcmp(d.timeUnit, 'years')
                    val = val([1 3]);
                else
                    val = val([2 3]);
                end
            end
            for k = 1:numel(d.lines)
                d.lines(k).XData = val{1}{k};
                d.lines(k).YData = val{2}{k};
            end
            % set axis props based on selection
            set(d.lines, 'Visible', 'on')
            switch sel
                case 'F-Phi'
                    d.Axes.XLabel.String = 'Storage capacity';
                    d.Axes.YLabel.String = 'Flow capacity';
                    d.Axes.XLim = [0 1];
                    d.Axes.YLim = [0 1];
                case 'Sweep vs time'
                    d.Axes.XLabel.String = ['Time [', d.timeUnit, ']'];
                    d.Axes.YLabel.String = 'Sweep';
                    d.Axes.YLim = [0 1];
                    d.Axes.XLim = [0 d.maxTOF];
                case 'Recovery factor vs Time'
                    d.Axes.XLabel.String = ['Time [', d.timeUnit, ']'];
                    d.Axes.YLabel.String = 'Recovery factor';
                    d.Axes.YLim = [0 1];
                    d.Axes.XLim = [0 d.maxTOF];
                case 'Residence time distribution'
                    d.Axes.XLabel.String = ['Time [', d.timeUnit, ']'];
                    d.Axes.YLabel.String = ['Recovery [1/', 'PVI',']'];
                    d.Axes.YLimMode = 'auto';
                    d.Axes.XLim = [0 d.maxTOF];
                case 'Oil Rates'
                    d.Axes.XLabel.String = ['Time [', d.timeUnit, ']'];
                    d.Axes.YLabel.String = 'Oil rate [m^3/day]';
                    d.Axes.YLimMode = 'auto';
                    d.Axes.XLim = [0 d.maxTOF];
                case 'Water Rates'
                    d.Axes.XLabel.String = ['Time [', d.timeUnit, ']'];
                    d.Axes.YLabel.String = 'Water rate [m^3/day]';
                    d.Axes.YLimMode = 'auto';
                    d.Axes.XLim = [0 d.maxTOF];
                case 'WCUT vs time'
                    d.Axes.XLabel.String = ['Time [', d.timeUnit, ']'];
                    d.Axes.YLabel.String = 'Water cut';
                    d.Axes.YLim = [0 1];
                    d.Axes.XLim = [0 d.maxTOF];
                case 'NPV vs time'
                    d.Axes.XLabel.String = ['Time [', d.timeUnit, ']'];
                    d.Axes.YLabel.String = 'NPV';
                    d.Axes.YLimMode = 'auto';
                    d.Axes.XLim = [0 d.maxTOF];
                otherwise
                    set(d.Axes, {'XLimMode', 'YLimMode'}, {'auto', 'auto'});
                    d.Axes.XLim = [0 d.maxTOF];
            end
        end
                
        
        function histPlotCallback(d, src, event)
            % handle histogram
            
            if ~strcmp(d.currentPlotType, 'hist')
                d.currentPlotType = 'hist';
                d.Menu.items{3}.ForegroundColor = d.textColor.inactive;
                d.Menu.items{4}.ForegroundColor = d.textColor.inactive;
                d.Menu.items{5}.ForegroundColor = d.textColor.inactive;
                d.Menu.items{6}.ForegroundColor = d.textColor.inactive;
                d.Menu.items{6}.items{1}.ForegroundColor  = d.textColor.inactive;
                d.Menu.items{6}.items{2}.ForegroundColor  = d.textColor.inactive;
                
                d.Menu.items{7}.ForegroundColor = d.textColor.active;
                set(d.markers, 'Visible', 'off')
                set(d.labels,  'Visible', 'off')
                set(d.lines,   'Visible', 'off')
                set(d.hists,   'Visible', 'off');
                d.memberSelectCallback();
            end
            
            set(d.hists, 'Visible', 'off');
            sel = d.selectionHistPlot;
            
            switch sel.selection
                case 'Porosity'
                    d.Axes.XLabel.String = 'Porosity';
                    d.Axes.YLabel.String = 'frequency';
                    unit = 1;  %dimentionles [0,1]
                case 'Perm.x'
                    d.Axes.XLabel.String = 'Perm.x [mD] ';
                    d.Axes.YLabel.String = 'frequency';
                    unit = (1e-3*darcy); %millidarcys
                case 'Perm.y'
                    d.Axes.XLabel.String = 'Perm.y [mD]';
                    d.Axes.YLabel.String = 'frequency';
                    unit = (1e-3*darcy); %millidarcys
                case 'Perm.z'
                    d.Axes.XLabel.String = 'Perm.z [mD]';
                    d.Axes.YLabel.String = 'frequency';
                    unit = (1e-3*darcy); %millidarcys
                otherwise
                    set(d.Axes, {'XLimMode', 'YLimMode'}, {'auto', 'auto'});
            end
            
            val = d.getDiagnosticsValues(sel.selection);
            
            hist_min  = d.Menu.items{7}.minValue*unit; %Transform to SI units;
            hist_max  = d.Menu.items{7}.maxValue*unit; %Transform to SI units;
            hist_bins = d.Menu.items{7}.binsValue;
            
            if ( d.Menu.items{7}.logSwitch == true)
                for k = 1:numel(d.hists)
                    cum_n_counts_log = val{k,3};
                    edges_log        = val{k,4};
                    
                    if any(edges_log)<= 0 % edges_log Data has zero and negatives values
                        warning('Logaritmic histogram is not available. There are zero values in data.');
                        d.hists(k).XData = edges_log/unit; % Transform to millidardys, in case of permeabilities
                        d.hists(k).YData = edges_log*0;
                    else % if there is non zero or negative we proceed with the logaricmic histogram
                        mini_aux =hist_min;
                        if (hist_min<=0)
                            mini_aux = hist_max*10^(-5);
                        end

                        EDGES_log        = linspace(log10(mini_aux),log10(hist_max),hist_bins+1); % hist_bins+1 =n_points
                        cum_HIST_N_log   = interp1(log10(edges_log),[0;cum_n_counts_log],EDGES_log') ;
                        N_log            = diff(cum_HIST_N_log);

                        [x_log,y_log]    = makeLinesFromHistogram(N_log,10.^EDGES_log);

                        d.hists(k).XData = x_log/unit; % Transform to millidardys, in case of permeabilities
                        d.hists(k).YData = y_log;
                    end
                end
            else
                for k = 1:numel(d.hists)
                    cum_n_counts = val{k,1};
                    edges        = val{k,2};
                    
                    EDGES      = linspace(hist_min,hist_max,hist_bins+1);
                    cum_HIST_N = interp1(edges,[0;cum_n_counts],EDGES') ;
                    N           = diff(cum_HIST_N);
                    [x,y] = makeLinesFromHistogram(N,EDGES);
                    
                    d.hists(k).XData = x/unit; % Transform to millidardys, in case of permeabilities   ;
                    d.hists(k).YData = y;
                end
            end
            
            % set axis props based on selection
            switch  d.Menu.items{7}.logSwitch
                case true
                    set(d.Axes, {'XLimMode', 'YLimMode'}, {'manual', 'auto'});
                    d.Axes.XLim = [d.Menu.items{7}.minValue d.Menu.items{7}.maxValue];
                    d.Axes.XScale = 'log';
                otherwise
                    set(d.Axes, {'XLimMode', 'YLimMode'}, {'manual', 'auto'});
                    d.Axes.XScale = 'linear';
                    d.Axes.XLim = [d.Menu.items{7}.minValue d.Menu.items{7}.maxValue];
            end
            set(d.hists, 'Visible', 'on')
        end
        
        
        
        % -----------------------------------------------------------------
        function modelViewerButtonCallback(d, src, event)
            launchDiagnosticsViewerForSelected(d);
        end
        % -----------------------------------------------------------------
        function thresholdCallback(d, src, event)
            if strcmp(d.currentPlotType, 'cross')
                d.crossPlotCallback(src, event);
            else
                d.linePlotCallback(src, event);
            end
        end
        % -----------------------------------------------------------------
        function fluidCallback(d, src, event)
            d.computePhaseRates();
            d.emptyCurrentDiagnostics();
            if strcmp(d.currentPlotType, 'cross')
                d.crossPlotCallback(src, event);
            else
                d.linePlotCallback(src, event);
            end
        end
        % -----------------------------------------------------------------
        function npvCallback(d, src, event)
            d.emptyCurrentDiagnostics({'NPV', 'NPVVsTime', 'maxNPV'})
            if strcmp(d.currentPlotType, 'cross')
                d.crossPlotCallback(src, event);
            else
                d.linePlotCallback(src, event);
            end
        end
        % -----------------------------------------------------------------
        function launchDiagnosticsViewerForSelected(d)
            nSel = nnz(d.memberSelection);
            if any(nSel)
                ii = find(d.memberSelection);
                d.diagnosticsViewerHandle = ...
                    DiagnosticsViewer(d.m, ii, 'modelNames',{d.names{ii}});
            end
        end
        % -----------------------------------------------------------------
        function regionSelectCallback(d, src, event)
            d.emptyCurrentDiagnostics();
            if strcmp(d.currentPlotType, 'cross')
                d.crossPlotCallback();
            else
                d.linePlotCallback();
            end
            
        end
        %------------------------------------------------------------------
        function d = computePhaseRates(d)
            showWaitbar = true;
            nMem = numel(d.validMembers);
            if showWaitbar
                wb = waitbar(0, 'Updating production rates for current fluid properties');
            end
            if ~d.useOriginalFluid
                fp = d.fluidProps;
                fluid = initSimpleADIFluid('phases', 'WO', 'mu', [fp.muw, fp.muo]*centi*poise,  ...
                    'smin', [fp.swc, fp.soc], 'n', fp.n);
            elseif isempty(d.Data.fluid)
                warning('Original fluid unavailable, rerun diagnostics to account for recent updates');
                close(wb);
                return;
            end
            
            for k = 1:nMem
                if showWaitbar
                    waitbar(k/nMem, wb);
                end
                if d.useOriginalFluid
                    fluid = d.Data.fluid{k};
                end
                d.Data.phaseRates{k} = computePhaseRatesFromRTD(d.Data.RTD{k}, fluid);
            end
            if showWaitbar, close(wb); end
        end
        
        %% main function for all computations ----------------------------
        function vals = getDiagnosticsValues(d, fld, range)
            if nargin < 3
                range = true;
            end
            if isempty(d.injectorIx) || isempty(d.producerIx)
                vals = cell(1, numel(d.validMembers));
                return;
            end
            
            memIx = reshape(d.validMembers, 1, []);
            nMem  = numel(memIx);
            switch matlab.lang.makeValidName(fld)
                case 'RTD'
                    if isempty(d.currentDiagnostics.RTD{1})
                        for k = 1:nMem
                            d.currentDiagnostics.RTD{k} = extractRegionRTD(d.Data.RTD{k}, d.injectorIx, d.producerIx);
                        end
                    end
                    vals = d.currentDiagnostics.RTD;
                case 'rTime'  % time vectors
                    if ~isempty(d.currentDiagnostics.rTime{1})
                        vals = d.currentDiagnostics.rTime;
                    else
                        r = d.getDiagnosticsValues('phaseRates');
                        vals = {cellfun(@(x)x.t/year,  r,   'UniformOutput', false), ...
                            cellfun(@(x)x.pvi,  r,   'UniformOutput', false)};
                        d.currentDiagnostics.rTime = vals;
                    end
                case 'ResidenceTimeDistribution'
                    if isempty(d.currentDiagnostics.ResidenceTimeDistribution{1})
                        rtd = d.getDiagnosticsValues('RTD');
                        d.currentDiagnostics.ResidenceTimeDistribution = ...
                            cellfun(@(x)x.values*x.volumes/x.allocations, rtd,   'UniformOutput', false);
                    end
                    vals = [d.getDiagnosticsValues('rTime'), {d.currentDiagnostics.ResidenceTimeDistribution}];
                case 'volume'
                    rtd = d.getDiagnosticsValues('RTD');
                    vals = cellfun(@(x)x.volumes, rtd);
                case 'allocation'
                    rtd = d.getDiagnosticsValues('RTD');
                    vals = cellfun(@(x)x.allocations, rtd);
                case 'WP'
                    if isempty(d.currentDiagnostics.WP{1})
                        for k = 1:nMem
                            d.currentDiagnostics.WP{k} = extractRegionWP(d.Data.WP{k}, d.injectorIx, d.producerIx);
                        end
                    end
                    vals = d.currentDiagnostics.WP;
                case 'F'
                    rtd = d.getDiagnosticsValues('RTD');
                    if isempty(d.currentDiagnostics.F{1})
                        for k = 1:nMem
                            if ~rtd{k}.allocations == 0
                                [d.currentDiagnostics.F{k}, d.currentDiagnostics.Phi{k}] = computeFandPhiFromDist(rtd{k});
                            else
                                [d.currentDiagnostics.F{k}, d.currentDiagnostics.Phi{k}] = deal([]);
                            end
                        end
                    end
                    if range
                        vals = getValuesWithinRange(d.currentDiagnostics.F, rtd, d.maxTOF, d.timeUnit, true);
                    else
                        vals = d.currentDiagnostics.F;
                    end
                case 'Phi'
                    rtd = d.getDiagnosticsValues('RTD');
                    if isempty(d.currentDiagnostics.Phi{1})
                        for k = 1:nMem
                            if ~rtd{k}.allocations == 0
                                [d.currentDiagnostics.F{k}, d.currentDiagnostics.Phi{k}] = computeFandPhiFromDist(rtd{k});
                            else
                                [d.currentDiagnostics.F{k}, d.currentDiagnostics.Phi{k}] = deal([]);
                            end
                        end
                    end
                    if range
                        vals = getValuesWithinRange(d.currentDiagnostics.Phi, rtd, d.maxTOF, d.timeUnit, true);
                    else
                        vals = d.currentDiagnostics.Phi;
                    end
                case 'F_Phi'
                    vals = {d.getDiagnosticsValues('Phi'), d.getDiagnosticsValues('F')};
                case 'SweepEfficiency'
                    if isempty(d.currentDiagnostics.SweepEfficiency{1})
                        F   = d.getDiagnosticsValues('F', false);
                        Phi = d.getDiagnosticsValues('Phi', false);
                        for k = 1:nMem
                            if ~isempty(F{k}) && ~isempty(Phi{k})
                                [Ev,tD] = computeSweep(F{k},Phi{k});
                            else
                                [Ev,tD] = deal(nan);
                            end
                            d.currentDiagnostics.SweepEfficiency{k} = struct('Ev', Ev(1:end-1), 'tD', tD(1:end-1));
                        end
                    end
                    vals = d.currentDiagnostics.SweepEfficiency;
                case 'SweepVsTime'
                    if isempty(d.currentDiagnostics.SweepVsTime{1})
                        rtd   = d.getDiagnosticsValues('RTD');
                        sweep = d.getDiagnosticsValues('SweepEfficiency');
                        d.currentDiagnostics.SweepVsTime = ...
                            {cellfun(@(x,y)x.tD*y.volumes/(y.allocations*year), sweep, rtd, 'UniformOutput', false), ...
                            cellfun(@(x)x.tD, sweep,   'UniformOutput', false), ...
                            cellfun(@(x)x.Ev, sweep,   'UniformOutput', false)};
                    end
                    vals = d.currentDiagnostics.SweepVsTime;
                case 'ModelNo'
                    vals = num2cell(memIx);
                case 'LorenzCoefficient'
                    F   = d.getDiagnosticsValues('F');
                    Phi = d.getDiagnosticsValues('Phi');
                    for k = 1:nMem
                        if ~isempty(F{k}) && ~isempty(Phi{k})
                            d.currentDiagnostics.LorenzCoefficient{k} = computeLorenz(F{k},Phi{k});
                        else
                            d.currentDiagnostics.LorenzCoefficient{k} = nan;
                        end
                    end
                    vals = d.currentDiagnostics.LorenzCoefficient;
                case 'Sweep'
                    sw   = d.getDiagnosticsValues('SweepEfficiency');
                    if strcmp(d.timeUnit, 'years')
                        rtd = d.getDiagnosticsValues('RTD');
                        tm  = cellfun(@(x,y)x.tD.*y.volumes/(y.allocations*year), sw, rtd, 'UniformOutput', false);
                    else
                        tm = cellfun(@(x)x.tD, sw, 'UniformOutput', false);
                    end
                    for k = 1:nMem
                        if ~isempty(sw{k}) && ~isempty(sw{k}.tD) && ~isempty(sw{k}.Ev)
                            d.currentDiagnostics.Sweep{k} = min(1, singleInterp(tm{k}, sw{k}.Ev, d.maxTOF));
                        else
                            d.currentDiagnostics.Sweep{k} = nan;
                        end
                    end
                    vals = d.currentDiagnostics.Sweep;
                case 'phaseRates'
                    if isempty(d.Data.phaseRates{1})
                        d.computePhaseRates();
                    end
                    if isempty(d.currentDiagnostics.WaterRates{1})
                        for k = 1:nMem
                            d.currentDiagnostics.phaseRates{k} = extractRegionRates(d.Data.phaseRates{k}, d.injectorIx, d.producerIx);
                        end
                    end
                    vals = d.currentDiagnostics.phaseRates;
                case 'WaterRates'
                    r = d.getDiagnosticsValues('phaseRates');
                    vals = [d.getDiagnosticsValues('rTime'), ...
                        {cellfun(@(x)x.ratesW*day, r, 'UniformOutput', false)}];
                case 'OilRates'
                    r = d.getDiagnosticsValues('phaseRates');
                    vals = [d.getDiagnosticsValues('rTime'), ...
                        {cellfun(@(x)x.ratesO*day, r, 'UniformOutput', false)}];
                case 'WCUTVsTime'
                    qw = d.getDiagnosticsValues('WaterRates');
                    qo = d.getDiagnosticsValues('OilRates');
                    vals = [qw(1:2), ...
                        {cellfun(@(w,o)w./(w+o), qw{3}, qo{3}, 'UniformOutput', false)}];
                case 'NPVVsTime'
                    if isempty(d.currentDiagnostics.NPVVsTime{1})
                        r  = d.getDiagnosticsValues('phaseRates');
                        op = d.npvProps;
                        dtvec = @(r)diff(r.t(:,1))./( (1+op.d/100).^(r.t(2:end,1)/year) );
                        npv = @(r)[0; cumsum(...
                            dtvec(r).*(op.ro*r.ratesO(2:end) - op.rwp*r.ratesW(2:end) ...
                            -op.rwi.*sum(r.ratesW(2:end) + r.ratesO(2:end), 2)))/stb];
                        d.currentDiagnostics.NPVVsTime = cellfun(npv, r, 'UniformOutput', false);
                    end
                    vals = [d.getDiagnosticsValues('rTime'), {d.currentDiagnostics.NPVVsTime}];
                case 'RecoveryVsTime'
                    if isempty(d.currentDiagnostics.RecoveryVsTime{1})
                        r  = d.getDiagnosticsValues('phaseRates');
                        d.currentDiagnostics.RecoveryVsTime = ...
                            cellfun(@(x)[0; cumsum(diff(x.t).*x.ratesO(2:end))], r, 'UniformOutput', false);
                    end
                    vals = [d.getDiagnosticsValues('rTime'), {d.currentDiagnostics.RecoveryVsTime}];
                case 'RecoveryFactorVsTime'
                    vals  = d.getDiagnosticsValues('RecoveryVsTime');
                    r     = d.getDiagnosticsValues('phaseRates');
                    vals{3} = cellfun(@(x,y)x/y.volumesO, vals{3}, r, 'UniformOutput', false);
                case 'RecoveryFactor'
                    v = d.getDiagnosticsValues('RecoveryFactorVsTime');
                    t = d.maxTOF;
                    if strcmp(d.timeUnit, 'years')
                        vals  = cellfun(@(x,y)singleInterp(x, y, t), v{1}, v{3}, 'UniformOutput', false);
                    else
                        vals  = cellfun(@(x,y)singleInterp(x, y, t), v{2}, v{3}, 'UniformOutput', false);
                    end
                case 'WCUT'
                    qo = d.getDiagnosticsValues('OilProduction');
                    qw = d.getDiagnosticsValues('WaterProduction');
                    vals = [qo(1:2), {cellfun(@(x,y)x./(x+y), qw{3}, qo{3}, 'UniformOutput', false)}];
                case 'NPVPerVolume'
                    v = d.getDiagnosticsValues('NPVVsTime');
                    a = num2cell(d.getDiagnosticsValues('volume'));
                    t = d.maxTOF;
                    if strcmp(d.timeUnit, 'years')
                        vals  = cellfun(@(x,y,z)singleInterp(x, y, t)/(z), v{1}, v{3}, a, 'UniformOutput', false);
                    else
                        vals  = cellfun(@(x,y,z)singleInterp(x, y, t)/(z), v{2}, v{3}, a, 'UniformOutput', false);
                    end
                case 'maxNPV'
                    if ~isempty(d.currentDiagnostics.maxNPV{1})
                        vals = d.currentDiagnostics.maxNPV;
                    else
                        v = d.getDiagnosticsValues('NPVVsTime');
                        vals = cellfun(@max, v{3}, 'UniformOutput', false);
                        d.currentDiagnostics.maxNPV = vals;
                    end
                case 'Porosity'
                    % Check if d.currentDiagnostics.Porosity is empty
                    if isempty(d.currentDiagnostics.Porosity{1,1})
                        for k = 1:nMem                            
                            eMem = d.validMembers(k);
                            d.currentDiagnostics.Porosity{k,1} = ...
                                d.m.diagnostics(1).other{eMem}.rock.poro.hist.n_counts;
                            d.currentDiagnostics.Porosity{k,2} = ...
                                d.m.diagnostics(1).other{eMem}.rock.poro.hist.edges;
                            d.currentDiagnostics.Porosity{k,3} = ...
                                d.m.diagnostics(1).other{eMem}.rock.poro.loghist.n_counts;
                            d.currentDiagnostics.Porosity{k,4} = ...
                                d.m.diagnostics(1).other{eMem}.rock.poro.loghist.edges;
                        end  
                    end
                    vals   = d.currentDiagnostics.Porosity;
                case 'Perm_x'
                    % Check if d.currentDiagnostics.Perm_x is empty or
                    % the number of bin has changed
                    if isempty(d.currentDiagnostics.Perm_x{1})
                        for k = 1:nMem
                            eMem = d.validMembers(k);
                            d.currentDiagnostics.Perm_x{k,1} = ...
                                d.m.diagnostics(1).other{eMem}.rock.perm.hist.n_counts(:,1);
                            d.currentDiagnostics.Perm_x{k,2} = ...
                                d.m.diagnostics(1).other{eMem}.rock.perm.hist.edges(:,1);
                            d.currentDiagnostics.Perm_x{k,3} = ...
                                d.m.diagnostics(1).other{eMem}.rock.perm.loghist.n_counts(:,1);
                            d.currentDiagnostics.Perm_x{k,4} = ...
                                d.m.diagnostics(1).other{eMem}.rock.perm.loghist.edges(:,1);
                         end
                    end
                    vals   = d.currentDiagnostics.Perm_x;
                case 'Perm_y'
                    % Check if d.currentDiagnostics.Perm_y is empty
                    if isempty(d.currentDiagnostics.Perm_y{1})
                        for k = 1:nMem
                            eMem = d.validMembers(k);
                            d.currentDiagnostics.Perm_y{k,1} = ...
                                d.m.diagnostics(1).other{eMem}.rock.perm.hist.n_counts(:,2);
                            d.currentDiagnostics.Perm_y{k,2} = ...
                                d.m.diagnostics(1).other{eMem}.rock.perm.hist.edges(:,2);
                            d.currentDiagnostics.Perm_y{k,3} = ...
                                d.m.diagnostics(1).other{eMem}.rock.perm.loghist.n_counts(:,2);
                            d.currentDiagnostics.Perm_y{k,4} = ...
                                d.m.diagnostics(1).other{eMem}.rock.perm.loghist.edges(:,2);
                        end
                    end
                    vals   = d.currentDiagnostics.Perm_y;
                case 'Perm_z'
                    % Check if d.currentDiagnostics.Perm_z is empty
                    if isempty(d.currentDiagnostics.Perm_z{1})
                        for k = 1:nMem
                            eMem = d.validMembers(k);
                            d.currentDiagnostics.Perm_z{k,1} = ...
                                d.m.diagnostics(1).other{eMem}.rock.perm.hist.n_counts(:,3);
                            d.currentDiagnostics.Perm_z{k,2} = ...
                                d.m.diagnostics(1).other{eMem}.rock.perm.hist.edges(:,3);
                            d.currentDiagnostics.Perm_z{k,3} = ...
                                d.m.diagnostics(1).other{eMem}.rock.perm.loghist.n_counts(:,3);
                            d.currentDiagnostics.Perm_z{k,4} = ...
                                d.m.diagnostics(1).other{eMem}.rock.perm.loghist.edges(:,3);
                        end
                    end
                    vals   = d.currentDiagnostics.Perm_z;
                otherwise
                    vals = cell(1, numel(d.validMembers));
                    warning('Unknown property: %s / %s', fld, matlab.lang.makeValidName(fld))
            end
            
        end
        
        
        %%-----------------------------------------------------------------
        function layout(d, ~, ~)
            [mw, sp] = deal(d.layoutParams.menuWidth, d.layoutParams.itemSpace);
            fip    = d.Figure.Position;
            mPos   = d.Menu.Position;
            mPos   = [5, fip(4)-mPos(4)-1, mw, mPos(4)];
            aPos   = [mw+2*sp, 2*sp, max(10, fip(3)-mw-3*sp), max(10, fip(4)-3*sp)];
            d.Menu.Position = mPos;
            d.Axes.Position = aPos;
        end
        %------------------------------------------------------------------
        function emptyCurrentDiagnostics(d, flds)
            if nargin < 2 % empty all
                flds = [{'WP', 'RTD', 'F', 'Phi', 'SweepEfficiency', 'recovery', 'rTime'}, ...
                    d.crossPlotSelections, d.linePlotSelections, d.histPlotSelections];
            end
            flds = cellfun(@matlab.lang.makeValidName, flds, 'UniformOutput', false);
            for k = 1:numel(flds)
                d.currentDiagnostics.(flds{k}) = cell(1, numel(d.validMembers));
            end
        end
        
        %------------------------------------------------------------------
        function [t, v] = collectValues(d, name, varargin)
            [t,v] = collectEnsembleGUIValues(d, name, varargin{:});
        end
    end
end
%--HELPER FUCNTIONS -------------------------------------------------------
%--------------------------------------------------------------------------
function vix = getValidMembers(m, ix, memberIx)
if isempty(memberIx)
    memberIx = 1:m.nMembers;
end
mix = false(m.nMembers,1);
mix(memberIx) = true;

vix = zeros(m.nMembers,1);
flds = {'WP', 'RTD', 'other'};
for k = 1:3
    vid = m.diagnostics(ix).(flds{k}).getValidIds();
    vix(vid) = vix(vid) + 1;
end
vix = find( (vix==3) & mix );

haveAll = numel(vix) == numel(memberIx);

if ~haveAll && ~isempty(vix)
    fprintf('Found computed diagnostics for %d of the requested %d realizations.\n', numel(vix), numel(memberIx));
end
end
%--------------------------------------------------------------------------
function data = loadAndProcessData(m, ix, mem)
data.WP  = m.diagnostics(ix).WP(mem);
data.RTD = m.diagnostics(ix).RTD(mem);
data.phaseRates  = cell(1, numel(mem));
other = m.diagnostics(ix).other(mem);
com = 0;
for k = 1:numel(mem)
    com = com + other{k}.wellCommunication/numel(mem);
end
data.wellCommunication = com;
data.fluid = [];
if isfield(other{1}, 'fluid')
    data.fluid = cellfun(@(x)x.fluid, other, 'UniformOutput', false);
end
end
%--------------------------------------------------------------------------
function rtdr = extractRegionRTD(rtd, iix, pix)
assert(strcmp(rtd.creator, 'computeRTD'), 'Need RTD from ''computeRTD'' for now')
rtdr = struct('pairIx', [1 1], 't', rtd.t(:,1), 'creator', rtd.creator);
if nargin == 1 || (isempty(iix)&&isempty(pix))
    [ix, iix, pix] = deal(':'); %sum all
else
    nw = max(rtd.pairIx);
    [liix, lpix] = deal(false(nw(1), 1), false(nw(2), 1));
    liix(iix) = true;
    lpix(pix) = true;
    ix = liix(rtd.pairIx(:,1)) & lpix(rtd.pairIx(:,2));
end

rtdr.volumes        = sum(rtd.volumes(ix));
rtdr.allocations    = sum(rtd.allocations(ix));
rtdr.injectorFlux   = sum(rtd.injectorFlux(iix));
rtdr.producerFlux   = sum(rtd.producerFlux(pix));
rtdr.pvi            = rtdr.t*rtdr.allocations/rtdr.volumes;

% each rtd is scaled such that integral equals alloc/injflux, hence
% multiply by corresponding injector-flux and scale back using region
% total injector flux
subRegInjFlux = rtd.injectorFlux(rtd.pairIx(ix,1));
rtdr.values   = (rtd.values(:, ix)*subRegInjFlux)/rtdr.injectorFlux;
end
%--------------------------------------------------------------------------
function rr = extractRegionRates(r, iix, pix)
assert(~isempty(r), 'This is a bug')
rr = struct('pairIx', [1 1], 't', r.t(:,1));
if nargin == 1 || (isempty(iix)&&isempty(pix))
    ix = ':'; %sum all
else
    nw = max(r.pairIx);
    [liix, lpix] = deal(false(nw(1), 1), false(nw(2), 1));
    liix(iix) = true;
    lpix(pix) = true;
    ix = liix(r.pairIx(:,1)) & lpix(r.pairIx(:,2));
end

rr.volumes        = sum(r.volumes(ix));
rr.volumesW       = sum(r.volumesW(ix));
rr.volumesO       = sum(r.volumesO(ix));
rr.allocations    = sum(r.allocations(ix));
rr.ratesW         = sum(r.ratesW(:, ix),2);
rr.ratesO         = sum(r.ratesO(:, ix),2);
if ~rr.volumes == 0
    rr.pvi = rr.t*rr.allocations/rr.volumes;
else
    rr.pvi = rr.t;
end
end
%--------------------------------------------------------------------------
function valsr = getValuesWithinRange(vals, rtd, maxTOF, timeUnit, normalize)
valsr = cell(1, numel(vals));
isPVI = strcmp(timeUnit, 'pvi');
for k = 1:numel(valsr)
    if ~isempty(vals{k})
        if isPVI
            ix = find(rtd{k}.pvi > maxTOF, 1, 'first');
        else
            ix = find(rtd{k}.t > maxTOF*year, 1, 'first');
        end
        if isempty(ix)
            ix = numel(rtd{k}.t);
        end
        valsr{k} = vals{k}(1:ix);
        if normalize
            valsr{k} = valsr{k}/valsr{k}(end);
        end
    end
end
end
% -----------------------------------------------------------------
% context menu helpers
function varargout = addAxesContextMenu(varargin)
nax = nargin;
varargout = cell(1,nax);
for k = 1:nax
    ax = varargin{k};
    assert(isa(ax, 'matlab.graphics.axis.Axes'));
    
    m = uicontextmenu('Parent', ax.Parent);
    uimenu(m, 'Label', 'Export to new figure', 'Callback', @(src, event)copyToFig(src, event, ax))
    
    ax.UIContextMenu = m;
    varargout{k} = ax;
end
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
% set(ax, 'FontSize', 12, 'Units', 'normalized');
% set(ax, 'OuterPosition', [0 0 1 1])
end
% -------------------------------------------------------------------------
function d = setupPlotOpts(d)
d.plotOpts.point.prop = {'LineStyle', 'Marker', 'MarkerSize',  'MarkerFaceColor', 'MarkerEdgeColor'};
d.plotOpts.point.val  = {'none',      'o',       12,           [.7 .7 .7],        [.2 .2 .2]};
d.plotOpts.pointSelect.prop   = {'MarkerFaceColor'};
d.plotOpts.pointSelect.val    = {'g'};
d.plotOpts.pointDeselect.prop = {'MarkerFaceColor'};
d.plotOpts.pointDeselect.val  = {[.7 .7 .7]};
d.plotOpts.line.prop = {'LineStyle', 'Marker', 'Color',    'lineWidth', 'Marker'};
d.plotOpts.line.val  = {'-',          '.'    , [.5 .5 .5],  .5,         'none'};
d.plotOpts.lineSelect.prop   = {'Color',    'lineWidth'};
d.plotOpts.lineSelect.val    = {'g', 4};
d.plotOpts.lineDeselect.prop = {'Color',    'lineWidth'};
d.plotOpts.lineDeselect.val  = {[.5 .5 .5],  1};
end
% -------------------------------------------------------------------------
function [xl, yl, corr] = fitLineSegment(x, y, xLim, yLim)
% fit line-segment to points x,y within box xLim, yLim
[xm, ym]   = deal(mean(x), mean(y));
[xb, yb]   = deal(x-xm, y-ym);
[xbxb, xbyb] = deal(xb'*xb, xb'*yb);
if xbxb > 0
    v = [xbxb, xbyb];
else
    v = [0, 1];
end
% any point p along line is of the form p = (xm, ym) + a*v
% get a-candidates
ac = [(xLim-xm)/v(1), (yLim-ym)/v(2)];
a  = [max(ac(ac<0)); min(ac(ac>0))];
[xl, yl] = deal(xm + a*v(1), ym + a*v(2));
if nargout > 1
    % compute correlation coefficent
    corr = abs(xbyb/(sqrt(xbxb)*norm(yb)));
end
end
%--------------------------------------------------------------------------
function v1 = singleInterp(t,v,t1)
% special-purpose interpolation for single input
ix = find(t1<t, 1, 'first');
if ~isempty(ix)
    if ix > 1
        a  = (v(ix)-v(ix-1))/(t(ix)-t(ix-1));
        v1 = v(ix)-a*(t(ix)-t1);
    else
        v1 = v(1);
    end
else
    v1 = v(end);
end
end
%--------------------------------------------------------------------------

function [x,y] = makeLinesFromHistogram(hist_N,hist_edges)
% This function create 2*n+2 points (x,y) to make a line from histogram
% variables of n bins (hist_N) and edges
j=1;
nn=length(hist_N);
x = zeros(1,nn*2+2); y=x;

x(j+0)=hist_edges(1); y(j+0)=0;
for i=1:nn
    x(j+1)=hist_edges(i); y(j+1)=hist_N(i);
    x(j+2)=hist_edges(i+1); y(j+2)=hist_N(i);
    j=j+2;
end
x(end)=hist_edges(nn+1); y(end)=0;
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
