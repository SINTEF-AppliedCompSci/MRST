classdef DiagnosticsViewer < handle
    properties
        Figure
        Axes3D
        Axes2DL
        Axes2DR
        colorBar
        colorHAx
        colormap3D = 'default';
        messageBar
        infoTextBox
        injColors
        prodColors
        WellPlot
        Data
        Menu
        Props
        Measures
        Allocation
        fdprops
        Distributions
        Patch
        maxTOF
        currentInteractionRegion = []
        currentFilterRegion      = []
        currentDiagnostics
        interactiveModes
        camlight
        outlineGrid
        info = []
        layoutParams = struct('menuWidth',        300, ...
            'distAxisHeight',   150, ...
            'itemSpace',         25, ...
            'includeToolTips', true);
    end
    
    methods
        function d = DiagnosticsViewer(arg1,arg2,varargin)
           
            % ------ Extract models ---------------------------------------
            [models,wells,states0] = extractInput(arg1, arg2);
            
            % ------ Load modules -----------------------------------------
            mrstModule add incomp mrst-gui
            
            % ------ Set options ------------------------------------------
            opt = struct('style',          'default',...
                         'modelNames',     [], ...
                         'maxTOF',         500*year, ...
                         'includeAverage', true);
            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            % ------ Create figure window ---------------------------------
            screensize = get( groot, 'Screensize' );
            wsize = .75*screensize(3:4);
            wsize = [screensize(3:4)-wsize*1.1, wsize];
            d.Figure = limitedToolbarFigure('Position', wsize);
            
            
            % ------ Setup messageBar --- ---------------------------------
            d.infoTextBox = uicontrol('Parent', d.Figure, ...
                'Style', 'text', ...
                'Position', [0 0 0 0], ...
                'Units', 'pixels', ...
                'String', '', ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment','left');
            %                 'BackgroundColor', [.5 .5 .5],...
            %                 'ForegroundColor',[1 1 1],...
            
            
            % ------ Setup messageBar --- ---------------------------------
            d.messageBar = uicontrol('Parent', d.Figure, ...
                'Style', 'text', ...
                'Position', [0 0 0 0], ...
                'Units', 'pixels', ...
                'String', '', ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment','left');
            %                 'BackgroundColor', [.5 .5 .5],...
            %                 'ForegroundColor',[1 1 1],...


            % ------ Setup model names --- ---------------------------------
            
            if isempty(opt.modelNames)
                modelNames = cell(1,numel(models));
                for i = 1:numel(models)
                    modelNames{i} = strcat('m',string(i));
                end
                opt.modelNames = modelNames;
            end
            
            if numel(opt.modelNames)>numel(models)
                opt.modelNames = opt.modelNames(1:numel(models));
            elseif numel(opt.modelNames) < numel(models)
                
                for i = 1:(numel(models)-numel(opt.modelNames))
                    modelNames{i} = strcat('m',string(i+numel(opt.modelNames)));
                end
                opt.modelNames = [opt.modelNames modelNames];
                d.messageBar.String = 'Fewer model names than models. Using generated names';
            end
            
            opt.modelNames = cellfun(@(x) char(x), opt.modelNames,'UniformOutput',false);
            
            
            % ------ Run diagnostics and get data -------------------------
            fprintf(1,'---------------------------------------\n') 
            fprintf(1,'Starting to compute diagnostics\n');
            d.Data = getDiagnosticsViewerData(models,wells,states0,extraOpt{:});
            d.maxTOF = opt.maxTOF;
            d.injColors = d.Data{1}.injColors;
            d.prodColors = d.Data{1}.prodColors;
            fprintf(1,'Diagnostics computed\n');
            
            
            % ------ Selectors for time steps and wells -------------------
            itemOpts = {'Parent', d.Figure, 'Visible','off', 'style', opt.style};
            
            selector3D.modsel = ...
                InputModelSelector('modelNames', opt.modelNames, itemOpts{:});
            m = getSelectedModelIndex(d,selector3D);
            
            % This assumes that the wells are the same in all model
            % realisations.
            selector3D.wsel = WellSelector(...
                'injectors', {d.Data{1}.diagnostics.WP.inj.name}, ...
                'producers', {d.Data{1}.diagnostics.WP.prod.name}, itemOpts{:});
            
            
            % ------ Selector for properties in 3D plot -------------------
            d.Props = struct(...
                'static',  struct('name', {{d.Data{m}.static.name}}, ...
                                  'limits', {{d.Data{m}.static.limits}}), ...
                'dynamic', struct('name', {{d.Data{m}.dynamic.name}}, ...
                                  'limits', {{d.Data{m}.dynamic.limits}}), ...
                'diagnostics', struct('name', {{...
                                  'TOF forward', 'TOF backward', ...
                                  'Residence time',  'Tracer forward', ...
                                  'Tracer backward', 'Tracer product', ...
                                  'Sweep regions',   'Drainage regions', ...
                                  'First arrival forward', ...
                                  'First arrival backward'}}, ...
                'limits', {{[0 1.001*d.maxTOF/year], ...
                            [0 1.001*d.maxTOF/year], ...
                            [0 2.001*d.maxTOF/year], [0 1], [0 1], [0 1], [], [], ...
                            [0 1.001*d.maxTOF/year], [0 1.001*d.maxTOF/year]}}), ...
                'computed', struct('name', {{}}, 'limits', {{}}));
            d.currentDiagnostics = emptyDiagnostics(d);
            selector3D.psel = PropertyDisplaySelector('props', d.Props, ...
                itemOpts{:}, 'includeLogSwitch', true, 'statisticsForAll',true);
            
            % ------ Selector for property filter in 3D -------------------
            selector3D.fsel = PropertyDisplaySelector( ...
                'Title', 'Filter display', 'includeFilter', true, ...
                'props', d.Props, 'includePlayer', false, ...
                itemOpts{:},...
                'includeEnableSwitch', true, 'includeLogSwitch', true, ...
                'startPlayCallback', @d.startPlayCallback, ...
                'stopPlayCallback', @d.stopPlayCallback);
            selector3D.fsel.Min = d.Data{m}.static(1).limits(1);
            selector3D.fsel.Max = d.Data{m}.static(1).limits(2);
            
            % ------ Selector for heterogeneity measures ------------------
            d.Measures = {{'none', 'F-Phi plot', 'Fractional recovery', ...
                           'Sweep efficiency', 'Lorenz coefficient'}};
            selector2D.msel = DynamicMeasureSelector('props', d.Measures, itemOpts{:});
            
            % ------ Selector for well allocation -------------------------
            d.Allocation = {{'none','Well connections', 'Injector volumes', ...
                             'Injector allocation', 'Injector profile', ...
                             'Producer volumes', 'Producer allocation', ...
                             'Producer profile', 'Well connection %'}};
            selector2D.asel = DynamicMeasureSelector(...
                'Title','Well allocations', ...
                'props', d.Allocation, ...
                'includeAvgSwitch', opt.includeAverage, itemOpts{:});
            
            % ------ Selector for summary output  -------------------------
            selector2D.ssel = WellSolSelector(d.Data{m}.wsdata.wellNames, ...
                d.Data{m}.wsdata.props, d.Data{m}.wsdata.times, itemOpts{:}, ...
                'showPlotWSButton',false, 'Title', 'Plot simulation output'); 
            
            % ------ Selector for RTD distribution  -----------------------
            selector2D.dsel = TracerSelector(itemOpts{:});
            
            % ------ Selector for fluid distribution ----------------------
            a = struct('name', {{d.Data{1}.dynamic.name}});
            b = struct('name',{a.name(2:end)});
            fdprops = cellfun(@(x) strcat({'Producer: '}, x(2:end),' volumes'),...
                              b.name, 'UniformOutput', false);
            d.fdprops = [{'none'}, {'Producer: all phases'}, horzcat(fdprops{:})];
            fdprops = cellfun(@(x) strcat({'Injector: '}, x(2:end),' volumes'),...
                              b.name, 'UniformOutput', false);
            d.fdprops = [d.fdprops, {'Injector: all phases'}, horzcat(fdprops{:})];
            selector2D.fdsel = ...
                TracerSelector('Title','Fluid Distributions',...
                               'props',d.fdprops,itemOpts{:});
            
            % ------ Create menu(s) ---------------------------------------
            % sub-menu for model selection
            m0 = UIMenu('Title', 'Model Selection', 'Parent', d.Figure, ...
                itemOpts{:}, 'items', {selector3D.modsel}, ...
                'collapseDirection', 'down');
            
            % sub-menu for property display (3d and 2d)
            m1 = UIMenu('Title', 'Property display selection', ...
                'Parent', d.Figure, itemOpts{:}, ...
                'items', {selector3D.psel, selector2D.ssel}, ...
                'collapseDirection', 'down');
            % sub-menu for region selection
            m2 = UIMenu('Title', 'Region selection', ...
                'Parent', d.Figure, itemOpts{:}, ...
                'items', {selector3D.wsel, selector3D.fsel}, ...
                'collapseDirection', 'down');
            % sub-menu for diagnostics plots
            m3 = UIMenu('Title', 'Diagnostics Analysis', ...
                'Parent', d.Figure, itemOpts{:}, ...
                'items', {selector2D.msel, selector2D.asel, selector2D.dsel}, ...
                'collapseDirection', 'down');
            % sub-menu for diagnostics plots
            m4 = UIMenu('Title', 'Multiphase Analysis', ...
                'Parent', d.Figure, itemOpts{:}, ...
                'items', {selector2D.fdsel}, ...
                'collapseDirection', 'down');
            
            % main menu
            d.Menu = UIMenu('Title', 'Menu', ...
                'Parent', d.Figure, itemOpts{:}, ...
                'items', {m0, m2, m1, m3, m4});
            
            
            % set different background color of each sub-menu: to this end,
            % select some light colors
            tt = [251 180 76; 58 152 216; 42 187 155; 252 121 122]./255;
            tt = [tt; tatarizeMap];
            for k =1:numel(d.Menu.items)
                d.Menu.items{k}.BackgroundColor = tt(k,:);
                try
                    d.Menu.items{k}.titleColor = 'w';
                catch
                    d.Menu.items{k}.ForegroundColor = 'w';
                end
            end
            
            % ------ Setup axes -------------------------------------------
            d.Axes3D  = axes('Parent', d.Figure, 'Units', 'pixels', 'CLimMode', 'manual');
            d.Axes2DL = axes('Parent', d.Figure, 'Units', 'pixels');
            d.Axes2DR = axes('Parent', d.Figure, 'Units', 'pixels');
            [d.Axes2DL, d.Axes2DR] = addAxesContextMenu(d.Axes2DL, d.Axes2DR);
            
            
            % ----- Collapse all menus except model selection menu --------
            m1.collapse = 1;
            m2.collapse = 1;
            m3.collapse = 1;
            m4.collapse = 1;
            for  f = reshape(fieldnames(selector2D),1, [])
                selector2D.(f{1}).collapse = 1;
            end
            
            % ------ Set figure callback functions ------------------------
            d.Figure.SizeChangedFcn        = @d.layout;
            d.Figure.WindowButtonDownFcn   = @d.resize;
            d.Figure.WindowButtonMotionFcn = @d.setMousePointer;
            d.Figure.WindowScrollWheelFcn  = @d.scrollMenu;
            
            % set zoom/pan/rotate
            d.interactiveModes           = setInteractiveModes(d.Axes3D);
            
            % ------ Show model with static property ----------------------
            d.Patch = CellDataPatch(d.Data{1}.G, d.Data{1}.static(1).values, ...
                                    'Parent', d.Axes3D, 'EdgeColor', [.3 .3 .3], ...
                                    'EdgeAlpha', .5, 'BackFaceLighting', 'lit');
            d.Patch = addPatchContextMenu(d.Patch);
            d.Figure.CurrentAxes = d.Axes3D;
            d.outlineGrid = plotGrid(d.Data{m}.G, 'FaceColor', 'none', ...
                                     'EdgeAlpha', 0.15, 'EdgeColor', [.3 .3 .3]);
            axis(d.Axes3D, 'tight', 'vis3d', 'off');
            d.Axes3D.ZDir = 'reverse';
            view(d.Axes3D, 3);
            daspect(d.Axes3D, [1 1 .2]);
 
            % ------ Add colorbar with histogram --------------------------
            d.colorBar = colorbar(d.Axes3D, 'WestOutside', 'AxisLocation', 'in');
            set(d.colorBar, 'Position',[.24 .8 .01 .19], 'Units', 'pixels', 'Xaxis', 'left');
            d.colorHAx = axes('Position',[.255 .8 .03 .19], 'Units', 'pixels');
            d.updateColorHist();
            
            % ------ Construct 2D plot axes -------------------------------
            axis(d.Axes2DL,'off');
            axis(d.Axes2DR,'off');
            
            % ------ Switch back to 3D plot axes and plot wells -----------
            % add extra toolbar stuff
            d.Figure.CurrentAxes = d.Axes3D;
            d.camlight    = light('Position', [1 1 -1], 'Style', 'infinite');
            d.camlight.Visible = 'off';
            d = addExtraTools(d);
            
            % finally wells
            d.WellPlot = WellPlotHandle(d.Data{m}.G, cellfun(@(data)data.wells, d.Data, 'UniformOutput', false), ...
                'Visible', 'off', 'Parent', d.Axes3D);
            for i=1:numel(d.WellPlot.producers)
                d.WellPlot.producers(i).label.FontSize = 8;
                d.WellPlot.producers(i).label.BackgroundColor = [.7 .7 .7];
            end
            for i=1:numel(d.WellPlot.injectors)
                d.WellPlot.injectors(i).label.FontSize = 8;
                d.WellPlot.injectors(i).label.BackgroundColor = [.7 .7 .7];
            end
            
            % ------ Set callbacks for 3D axes ----------------------------
            selector3D.modsel.Callback = @(src, event) ...
                d.modelCallback(src, event, selector3D, selector2D);
            selector3D.modsel.modelViewerButtonCallback = @(src, event) ...
                d.modelViewerButtonCallback(src, event, selector3D);
            selector3D.psel.Callback = @(src, event) ...
                d.displayPropCallback(src, event, selector3D);
            selector3D.wsel.Callback = @(src, event) ...
                d.interactionRegionCallback(src, event, selector2D, selector3D);
            selector3D.fsel.Callback = @(src, event) ...
                d.filterPropCallback(src, event, selector3D);
            
            d.infoTextBox.Callback = @(src,event) ...
                d.modelInfoDisplayCallback(src, event, selector3D);
            
            % ------ Set callbacks for 2D axes ----------------------------
            selector2D.msel.Callback = @(src, event) ...
                d.measureCallback(src, event, selector2D, selector3D);
            selector2D.asel.Callback = @(src, event) ...
                d.allocationCallback(src, event, selector2D, selector3D);
            
            selector2D.ssel.Callback = @(src, event) ...
                d.wellSolCallback(src, event, selector2D, selector3D);
            selector2D.ssel.plotWellSolCallback = @(src, event) ...
                d.plotWellSolCallback(src, event, selector3D);
            
            selector2D.ssel.regCallback = @(src, event) ...
                d.selectWellsForSummary(src, event, selector2D, selector3D);
            selector2D.dsel.Callback    = @(src, event) ...
                d.distributionCallback(src, event, selector2D, selector3D);
            selector2D.fdsel.Callback   = @(src, event) ...
                d.fluidDistributionCallback(src, event, selector2D, selector3D);
            
            % ------ Do initial callback for 3D axes ----------------------
            selector3D.psel.Callback([], [])
            d.layout();
            
            d.infoTextBox.Callback([],[])
        end
        
        %% ---------- MAIN CALLBACKS --------------------------------------
        function modelCallback(d, src, event, s3, s2)
            m = s3.modsel.ix;
            d.currentDiagnostics = emptyDiagnostics(d);
            if isempty(m) % disable all controls except
                s3.modsel.Enable = 'on';
                s3.wsel.Enable = 'off';   s3.wsel.collapse = 1;
                s2.msel.Enable = 'on';    s2.msel.collapse = 1;
                s2.asel.Enable = 'off';   s2.asel.collapse = 1;
                s2.dsel.Enable = 'off';   s2.dsel.collapse = 1;
                set([s3.psel.typePopup, s3.psel.propPopup], 'Enable', 'on');
                s3.psel.typeIx = 1;
                s3.psel.propIx = 1;
                s3.psel.propPopup.String = {d.Data{1}.static.name};
                d.displayPropCallback(src, event, s3);
                d.modelInfoDisplayCallback(src, event, s3);
            else
                d.messageBar.String = '';
                s3.wsel.Enable  = 'on';  s3.wsel.collapse = 0;
                s2.msel.Enable  = 'on';  s2.msel.collapse = 0;
                s2.dsel.Enable  = 'on';
                s2.asel.Enable  = 'on';
                s2.rsel.Enable  = 'on';
                s2.fdsel.Enable  = 'on';
                s3.fsel.collapse = 1;
                s3.wsel.collapse = 0;
                if numel(m) == 1
                    % some stat-stuff
                end
                d.WellPlot.visibleCases = m;
                d.displayPropCallback(src, event, s3);
                d.interactionRegionCallback(src, event, s2, s3);
                s3.wsel.communicationMatrix = getSelectedModelCommunicationMatrix(d,s3);
                d.modelInfoDisplayCallback(src, event, s3);
            end
        end
        % -----------------------------------------------------------------
        function displayPropCallback(d, src, event, s3) %#ok<*INUSL>
            m      = getSelectedModelIndex(d,s3);
            propIx = s3.psel.propIx;
            numM   = numel(m);
            switch s3.psel.typeIx
                case 1  % static property
                    if numM == 1
                        displayVals = d.Data{m}.static(propIx).values;
                        lims = d.Data{m}.static(propIx).limits;
                    else
                        d.messageBar.String  = '';
                        stat = s3.psel.statPopup.String{s3.psel.statIx};
                        cval = zeros(numel(d.Data{m(1)}.static(propIx).values),numM);
                        for j = 1:numM
                            cval(:,j) = d.Data{m(j)}.static(propIx).values;
                        end
                        displayVals = computeStatistic(cval, stat);
                        lims = [min(displayVals) max(displayVals)];
                    end
                    
                case 2  % dynamic property
                    if numM == 1
                        displayVals = d.Data{m}.dynamic(propIx).values;
                        lims = d.Data{m}.dynamic(propIx).limits;
                    else % multiple models selected
                        d.messageBar.String  = '';
                        stat = s3.psel.statPopup.String{s3.psel.statIx};
                        cval = zeros(numel(d.Data{m(1)}.dynamic(propIx).values),numM);
                        for j = 1:numM
                            cval(:,j) = d.Data{m(j)}.dynamic(propIx).values;
                        end
                        displayVals = computeStatistic(cval, stat);
                        lims = [min(displayVals) max(displayVals)];
                    end
                    
                case 3   % diagnostics property
                    prop = s3.psel.propPopup.String{propIx};
                    
                    [d, cval, lims, flag] = ...
                        extractSelectedDiagnostics(d, prop, s3.modsel, s3.wsel);
                    if numM == 1 || flag
                        displayVals = cval;
                        s3.psel.statIx = 1;
                    else % multiple time steps
                        d.messageBar.String = '';
                        stat = s3.psel.statPopup.String{s3.psel.statIx};
                        displayVals = computeStatistic(cval, stat, prop); % include prop
                    end
                    
                case 4   % computed property
                    if ~isempty(d.Props.computed.name)
                        displayVals = d.Data{m(1)}.computed(propIx).values;
                        lims        = d.Data{m(1)}.computed(propIx).limits;
                    else % no computed properties return to static
                        d.messageBar.String = 'No computed values available';
                        s3.psel.typeIx = 1;
                        s3.psel.propPopup.String = {d.Data{1}.static(1).name};
                        return
                    end
            end
            isLog = false;
            if s3.psel.logSwitch
                [lims, displayVals, flag] = makeSafeForLog(lims, displayVals, 5);
                if ~flag % values not good for log-plot, reset switch
                    s3.psel.logSwitchBox.Value = 0;
                else
                    isLog = true;
                    [lims, displayVals] = deal(log10(lims), log10(displayVals));
                end
            end
            d.Patch.colorData = displayVals;
            if diff(lims) ~= 0
                d.Axes3D.CLim = lims;
            else
                d.Axes3D.CLimMode = 'auto';
                lims = d.Axes3D.CLim;
            end
            d.updateColorHist();
            d.updateColorBar(s3, lims, isLog);
            d.modelInfoDisplayCallback(src, event, s3);
        end
        % -----------------------------------------------------------------
        function filterPropCallback(d, src, event, s3)
            m      = getSelectedModelIndex(d,s3);
            propIx = s3.fsel.propIx;
            numM   = numel(m);
            if ~strcmp(s3.fsel.playMode, 'stop')
                d.Patch.value = s3.fsel.maxValue;
            else
                if strncmp(src.Style,'check',5) && ...
                        src.Value == 0 && strncmp(src.String,'Enable',6)
                    d.currentFilterRegion = [];
                    if ~isempty(d.currentInteractionRegion)
                        d.Patch.cells = d.currentInteractionRegion;
                    else
                        d.Patch.cells = 1:d.Data{1}.G.cells.num;
                    end
                    return
                end
                switch s3.fsel.typeIx
                    case 1  % static property
                        if numM>1, m = m(1); end
                        fval = d.Data{m}.static(propIx).values;
                        
                    case 2  % dynamic property
                        if numM == 1
                            fval = d.Data{m}.dynamic(propIx).values;
                        else  % multiple time steps
                            d.messageBar.String = '';
                            stat = s3.fsel.statPopup.String{s3.fsel.statIx};
                            for j = 1:numM
                                fval = zeros(numel(d.Data{m(j)}.dynamic(propIx).values),numM);
                                fval(:,j) = d.Data{m(j)}.dynamic(propIx).values;
                            end
                            fval = computeStatistic(fval, stat);
                        end
                        
                    case 3   % diagnostics property
                        prop = s3.fsel.propPopup.String{propIx};
                        [d, fval] = extractSelectedDiagnostics(d, prop, s3.modsel, s3.wsel);
                        if numM > 1 % multiple models steps
                            stat = s3.fsel.statPopup.String{s3.fsel.statIx};
                            fval = computeStatistic(fval, stat, prop); % include prop
                        end
                        
                    case 4 % computed prop
                        if ~isempty(d.Props.computed.name)
                            fval = d.Data{m}.computed(propIx).values;
                        else % no computed properties return to static
                            d.messageBar.String = 'No computed values available';
                            s3.fsel.typeIx = 1;
                            s3.fsel.propPopup.String = {d.Data{1}.static(1).name};
                            return
                        end
                end
                d.Patch.filterData = fval;
                d.currentFilterRegion = and(s3.fsel.minValue <= fval, s3.fsel.maxValue >= fval);
                if isempty(d.currentInteractionRegion)
                    d.Patch.cells = d.currentFilterRegion;
                else
                    d.Patch.cells =  and(d.currentFilterRegion, d.currentInteractionRegion);
                end
                d.updateColorHist();
            end
        end
        % -----------------------------------------------------------------
        function interactionRegionCallback(d, src, event, s2, s3)
            d.WellPlot.visibleInjectors = s3.wsel.injectorIx;
            d.WellPlot.visibleProducers = s3.wsel.producerIx;
            if isprop(src,'Style') && all(strncmp(src.Style,'checkbox',6))
                return;
            end
            % clear current diagnostics
            % this should only be done if well-selection has
            % changed, not if threshold is changed -> split into two
            % callbacks
            d.currentDiagnostics = emptyDiagnostics(d);
            [d, val] = extractSelectedDiagnostics(d, 'Tracer product', s3.modsel, s3.wsel);
            if numel(s3.modsel.ix) > 1 % do mean
                val = computeStatistic(val, 'mean');
            end
            d.currentInteractionRegion = val >= s3.wsel.threshold;
            if isempty(d.currentFilterRegion)
                d.Patch.cells = d.currentInteractionRegion;
            else
                if isprop(src,'Style') && all(strncmp(src.Style,'listbox',6))
                    d.filterPropCallback(src, event, s3);
                else
                    d.Patch.cells = ...
                        and(d.currentFilterRegion, d.currentInteractionRegion);
                end
            end
            if s3.psel.typeIx == 3 % reset display
                d.displayPropCallback(src, event, s3);
            end
            d.updateColorHist();
            % when interaction region has changed, set small axis color to
            % gray to indicate not up-to-date
            d.Axes2DL.Color = [.85 .85 .85];%d.Figure.Color;
            d.Axes2DR.Color = [.85 .85 .85];%d.Figure.Color;
            if s2.ssel.regionSwitch.Value == 1
                s2.ssel.regCallback(src, event)
            end
        end
        % -----------------------------------------------------------------
        function measureCallback(d, src, event, s2, s3)
            if s2.msel.panelNo==1
                ax = d.Axes2DL;
                d.resetSelectors(s2, 'msel','leftIx');
            else
                ax = d.Axes2DR;
                d.resetSelectors(s2, 'msel', 'rightIx');
            end
            cla(ax, 'reset');
            resetValue = false;
            switch src.Value
                case 1
                    axis(ax,'off');
                case 2
                    if isempty(s3.modsel.ix)
                        axis(ax,'off'); resetValue = true;
                    else
                        d.showFPhi(ax, s3.modsel, 1, s3.wsel);
                    end
                case 3
                    if isempty(s3.modsel.ix)
                        axis(ax,'off'); resetValue = true;
                    else
                        d.showFracRecovery(ax, s3.modsel, 1, s3.wsel);
                    end
                case 4
                    if isempty(s3.modsel.ix)
                        axis(ax,'off'); resetValue = true;
                    else
                        d.showSweep(ax, s3.modsel, 1, s3.wsel);
                    end
                case 5
                    if isempty(s3.modsel.ix)
                        axis(ax,'off'); resetValue = true;
                    else
                        d.showLorenz(ax, s3.modsel, 1, s3.wsel);
                    end
                otherwise
                    disp('functionality not implemented yet');
            end
            if resetValue
                if s2.msel.panelNo==1
                    s2.msel.leftIx = 1;
                else
                    s2.msel.rightIx = 1;
                end
            end
        end
        % -----------------------------------------------------------------
        function allocationCallback(d, src, event, s2, s3)
            if s2.asel.panelNo==1
                ax = d.Axes2DL;
                d.resetSelectors(s2, 'asel', 'leftIx');
            else
                ax = d.Axes2DR;
                d.resetSelectors(s2, 'asel', 'rightIx');
            end
            cla(ax, 'reset');
            showAllocation(d, ax, s2, s3)
        end
        % -----------------------------------------------------------------
        function distributionCallback(d, src, event, s2, s3)
            [dsel, modsel, wsel] = deal(s2.dsel, s3.modsel, s3.wsel);
            if isempty(modsel.ix) || dsel.extendTime <= 0, return, end
            if s2.dsel.panelNo==1
                [ax, ix] = deal(d.Axes2DL, dsel.leftIx);
            else
                [ax, ix] = deal(d.Axes2DR, dsel.rightIx);
            end
            cla(ax, 'reset');
            if numel(wsel.injectorIx) ~= 1 || numel(modsel.ix) ~= 1
                d.messageBar.String = 'Please select  one injector and  one model.';
                return
            else
                d.messageBar.String = '';
                m = modsel.ix;
                [iIx, pIx] = deal(s3.wsel.injectorIx, s3.wsel.producerIx);
                if isempty(pIx)
                    pIx = 1:d.WellPlot.nProd;%    1:numel(d.WellPlot.producers);
                end
                switch ix
                    case 1
                        axis(ax, 'off')
                        return
                    case 2 % estimate
                        PORVIx = arrayfun(@(x) strcmp(x.name,'PORV'),d.Data{m}.static);
                        PORV = d.Data{m}.static(PORVIx).values;
                        
                        dist = estimateRTD(PORV, d.Data{m}.diagnostics.D, ...
                            d.Data{m}.diagnostics.WP, ...
                            'injectorIx', iIx, 'producerIx', pIx);
                    case 3 % compute
                        PORVIx = arrayfun(@(x) strcmp(x.name,'PORV'),d.Data{m}.static);
                        PORV = d.Data{m}.static(PORVIx).values;
                        dist = computeRTD(d.Data{m}.states, d.Data{m}.Gs, PORV, ...
                            d.Data{m}.diagnostics.D, d.Data{m}.diagnostics.WP, ...
                            d.Data{m}.wells, ...
                            'injectorIx', iIx, 'producerIx', pIx);
                end
                for k = 1:numel(pIx)
                    line(ax, dist.t(:,k)/year, dist.values(:,k), ...
                        'LineWidth', 2, 'Color', d.Data{m}.prodColors(pIx(k),:));
                end
                ylabel(ax, 'Tracer rate');
                ax.XLabel.String = 'TOF distance in years';
                %wn = arrayfun(@(x)x.label.String, d.WellPlot.producers(pIx), 'UniformOutput', false);
                wn = d.WellPlot.producerNames(pIx);
                legend(ax, wn, 'Location','northeast', 'Interpreter', 'none')
                set(ax, 'FontSize', 10)
                set(ax, 'XLim', [0, s2.dsel.extendTime]);
            end
        end
        % -----------------------------------------------------------------
        function wellSolCallback(d, src, event, s2, s3, ax)           
            m = getSelectedModelIndex(d,s3);           
            if nargin < 6
                if s2.ssel.panelNo == 1
                    ax = d.Axes2DL;
                else
                    ax = d.Axes2DR;
                end
            end
            [nms, prps] = deal(s2.ssel.curNames, s2.ssel.curProps);
            
            if numel(prps)>1
                d.messageBar.String = strcat("Multiple selected properties not allowed. Plotting ",...
                    prps(1)," data.");
            end
            
            propIx = s2.ssel.propIx;
            leg={};
            nameIx = zeros(numel(nms),1);
            for l = 1:numel(nms)
                nameIx(l) = find(strcmp(d.Data{1}.wsdata.wellNames, nms(l)));
            end
            
            if isempty(prps)
                return;
            end
            if nargin < 6,  cla(ax, 'reset'); end
            
            hold(ax, 'on')
            values = zeros(numel(m),numel(nms));
            
            if numel(m)>1
                for k = numel(m):-1:1
                    for l = 1:numel(nms)
                        values(k,l) = d.Data{m(k)}.wsdata.(prps{1})(:,nameIx(l));
                    end
                    groupnames{k} = s3.modsel.modelNames{m(k)};
                end
                leg = {};
                for l = numel(nms):-1:1
                    leg(l) = {[d.Data{1}.wsdata.wellNames{nameIx(l)}]};
                end
            else
                for l = 1:numel(nms)
                    values(1,l) = d.Data{m}.wsdata.(prps{1})(:,nameIx(l));
                end
                groupnames= nms;
            end
            
            colors = [d.Data{1}.injColors; d.Data{1}.prodColors];
            
            h=bar(ax, categorical(groupnames),values,'grouped','FaceColor','flat');
            hold(ax,'on');
            if numel(m)>1
                for k = 1:size(values,2)
                    h(k).CData = colors(nameIx(k),:);
                end
            else
                h.CData = colors(nameIx,:);
            end
            if ~isempty(leg)
                 legend(ax, leg, 'Interpreter', 'none')
            end
            ylabel(ax, {[prps{1}, ' [', d.Data{1}.wsdata.units{propIx(1)}, ']']})
        end
        
        function modelViewerButtonCallback(d, src, event, s3)
            m = getSelectedModelIndex(d,s3);
            if numel(m) > 16
                d.messageBar.String = 'Please select fewer than 16 models.';
                return
            else
                d.messageBar.String = '';
            end
            
            % Get data
            cval = cell(1,numel(m));
            switch s3.psel.typeIx
                case 1
                    prop = 'static';
                    for j = 1:numel(m)
                        cval{j} = d.Data{m(j)}.(prop)(s3.psel.propIx).values;
                    end
                    plottitle = d.Data{1}.(prop)(s3.psel.propIx).name;
                case 2
                    prop = 'dynamic';
                    for j = 1:numel(m)
                        cval{j} = d.Data{m(j)}.(prop)(s3.psel.propIx).values;
                    end
                    plottitle = d.Data{1}.(prop)(s3.psel.propIx).name;
                case 3
                    prop = s3.psel.propPopup.String{s3.psel.propIx};
                    for j = 1:numel(m)
                        [d, vals , ~, flag] = ...
                            extractSelectedDiagnosticsSingleModel(d, prop, m(j), s3.wsel);
                        if flag
                            cval{j} = vals;
                        else
                            for i = 1:numel(m)
                                cval{i} = vals(:,i);
                            end
                            break
                        end
                    end
                    plottitle = prop;
                case 4
                    prop = 'computed';
                    plottitle = prop;
            end
            
            % Set colorbar limits
            if size(cval{1},2)==3
                clims = [0 1];
            else
                clims = [min(cellfun(@min,cval)) max(cellfun(@max,cval))];
            end
            
            if s3.psel.logSwitch
                isLog = true;
                for i = 1:numel(cval)
                    [clims, cval{i}, flag] = makeSafeForLog(clims, cval{i}, 5);
                    if ~flag % values not good for log-plot, reset switch
                        isLog = false;
                        return
                    end
                end
                if isLog
                    cval  = cellfun(@log10, cval, 'UniformOutput', false);
                    clims = log10(clims);
                end
            else
                isLog = false;
            end
            
            % Plot figure
            plotMultipleLinkedModels( d.Data{m(1)}.G, cval, clims,...
                d.Patch.cells, s3.modsel.modelNames(m), isLog, s3, d.Data{1}.wells,...
                d.injColors, d.prodColors, plottitle, d.Axes3D.View, d.Axes3D.CLim,...
                'camlight', d.camlight, 'outlineGrid', d.outlineGrid);
            % 'filterMinVal',s3.fsel.minValue,'filterMaxVal',s3.fsel.maxValue);   
        end
        
        function plotWellSolCallback(d, src, event, s3)
            m = getSelectedModelIndex(d,s3);
            plotWellSols({d.Data{m(1)}.wellSol},d.Data{m(1)}.wsdata.times);
        end
        % -----------------------------------------------------------------
        function selectWellsForSummary(d, src, event, s2, s3)
            prd = s3.wsel.prodSelector.String(s3.wsel.producerIx);
            inj = s3.wsel.injSelector.String(s3.wsel.injectorIx);
            s2.ssel.setWellSubset([prd(:);inj(:)]);
        end
        
        % -----------------------------------------------------------------
        function fluidDistributionCallback(d, src, event, s2, s3)           
            if s2.fdsel.panelNo==1
                [ax, ix] = deal(d.Axes2DL, s2.fdsel.leftIx);
            else
                [ax, ix] = deal(d.Axes2DR, s2.fdsel.rightIx);
            end
            cla(ax, 'reset');
            d.messageBar.String = "";
            
            % If first menu option ('none'): clear axis and return
            if ix==1
                axis(ax, 'off');
                return;
            end
            % If more than one model selected, clear axis and return
            if numel(s3.modsel.ix)>1
                d.messageBar.String = strcat("Multiple models selected. Please select only one model.");
                return
            end
            
            if ix<5
                wellIx = s3.wsel.producerIx;
                wtype = 'producer';
            else
                wellIx = s3.wsel.injectorIx;
                wtype  = 'injector';
            end
            
            if numel(wellIx)<1
                d.messageBar.String = ...
                    ['No ',wtype,' selected. Please select one ',wtype,'.'];
                return
            elseif numel(wellIx)>1
                d.messageBar.String = ...
                    ['Multiple ',wtype,'s selected. Please select one ',wtype,'.'];
                return
            else
                d.messageBar.String = "";
            end
            switch ix
                case {2,5}
                    % Plot fluid  arrivals
                    hold(ax,'off');
                    d.PlotTOFArrival(ax, s3.modsel, 1, wellIx, wtype, s2.fdsel.extendTime);
                case {3,4}
                    d.PlotPhaseTOFArrivals(ax, s3.modsel, 1, s3.wsel, ...
                                           wtype, ix-2, s2.fdsel.extendTime);
                case {6,7}
                    d.PlotPhaseTOFArrivals(ax, s3.modsel, 1, s3.wsel, ...
                                           wtype, ix-5, s2.fdsel.extendTime);
                otherwise
                    error('fluidDistributionCallback: something is wrong');
            end
            
        end
        % -----------------------------------------------------------------
        function modelInfoDisplayCallback(d, src, event, s3)     
            txt = d.getModelInfoString(s3);
            nlines = numel(strsplit(txt,'\n'));
            d.infoTextBox.Position(4) = 18*nlines;
            d.infoTextBox.String = txt;        
        end
        % -----------------------------------------------------------------
        function startPlayCallback(d, src)
            % set mouse pointer to wait-mode
            fsel = src.Parent.UserData;
            mtfnc  = d.Figure.WindowButtonMotionFcn;
            curpnt = d.Figure.Pointer;
            d.Figure.WindowButtonMotionFcn = '';
            d.Figure.Pointer = 'watch';
            pause(0.001); % change pointer now
            set(d.Axes3D, {'YLimMode', 'XLimMode', 'ZLimMode'}, {'manual', 'manual', 'manual'});
            d.Patch.logScale = fsel.logSwitch;
            d.Patch.playMode = true;
            % disable all exept playbar/slider
            cellfun(@(x)set(x, 'Enable', 'off'), d.Menu.items)
            set(fsel.playBar, 'Enable', 'on');
            fsel.maxSlider.Enable = 'on';
            %reset pointer
            d.Figure.Pointer = curpnt;
            d.Figure.WindowButtonMotionFcn = mtfnc;
        end
        % -----------------------------------------------------------------
        function stopPlayCallback(d, ~, ~)
            d.Patch.playMode = false;
            set(d.Axes3D, {'YLimMode', 'XLimMode', 'ZLimMode'}, {'auto', 'auto', 'auto'});
            % enable all
            cellfun(@(x)set(x, 'Enable', 'on'), d.Menu.items)
        end
        % -----------------------------------------------------------------
        %%  -------- LAYOUT CALLBACKS -------------------------------------
        function layout(d, ~, ~)
            [mw, ah, sp] = deal(d.layoutParams.menuWidth, d.layoutParams.distAxisHeight, d.layoutParams.itemSpace);
            fip    = d.Figure.Position;
            mPos   = d.Menu.Position;
            mPos   = [5, fip(4)-mPos(4)-1, mw, mPos(4)];
            mPos(3) = max(mPos(3), sp);
            aPos2D = [mw+sp, 2.5*sp, (fip(3)-mw-4*sp)/2, ah];
            aPos2D(4) = max(aPos2D(4), 0);
            aPos2DL = aPos2D; aPos2DL([1 3]) = aPos2DL([1 3]) + [sp -sp];
            aPos2DR = aPos2D; aPos2DR(1) = aPos2D(1)+aPos2D(3)+2*sp;
            aPos3D = [mw+sp, 2.5*sp+ah, fip(3)-mw-2*sp, fip(4)-3*sp-ah];
            cbh    = max(50, min(300, fip(4)-2*sp));
            cbw    = 37;

            % Colorbar left, next to menu
            aPosCB = [mw+2*sp,       fip(4)-cbh-sp, cbw, cbh];
            aPosHA = [mw+2*sp+cbw+5, fip(4)-cbh-sp, cbw, cbh];
            
            % Colorbar right
            d.Menu.Position     = mPos;
            d.messageBar.Position = [sum(mPos(1,3))+2, 1, fip(3)-sum(mPos(1,3))-2, 20];
            d.infoTextBox.Position = [5, 5, mPos(3), 200];
            d.Axes2DL.Position  = aPos2DL;
            d.Axes2DR.Position  = aPos2DR;
            d.Axes3D.Position   = max(0,aPos3D);
            d.colorBar.Position = aPosCB;
            d.colorHAx.Position = aPosHA;
            d.infoTextBox.Callback([],[]);
        end
        % -----------------------------------------------------------------
        function resize(d, src, ~)
            mw = d.layoutParams.menuWidth; 
            ah = d.layoutParams.distAxisHeight;
            sp = d.layoutParams.itemSpace;
            p  = src.CurrentPoint;
            if p(1) <= mw % inside menu, deal with it there
                % do nothing
            elseif p(1) < mw+sp   % drag menu width
                motionFcn = d.Figure.WindowButtonMotionFcn;
                offset    = p(1)-mw;
                src.WindowButtonMotionFcn = {@dragHorisontally, offset};
                src.WindowButtonUpFcn     = {@clearMotionFcn, motionFcn};
            elseif (p(2)>ah+2*sp) && (p(2)<ah+3*sp)
                motionFcn = d.Figure.WindowButtonMotionFcn;
                offset    = p(2)-ah-sp;
                src.WindowButtonMotionFcn = {@dragVertically, offset};
                src.WindowButtonUpFcn     = {@clearMotionFcn, motionFcn};
            end
            
            function dragHorisontally(src, ~, offset)
                p = src.CurrentPoint - [offset, 0];
                if abs(mw - p(1)) > 5
                    d.layoutParams.menuWidth = max(p(1), 50);
                    d.Menu.Position(3) = d.layoutParams.menuWidth;
                    %d.layout();
                end
            end
            
            function dragVertically(src, ~, offset)
                p = src.CurrentPoint - [0 offset];
                if abs(ah - p(2)) > 2
                    d.layoutParams.distAxisHeight = max(p(2)-sp, sp);
                    d.layout();
                end
            end
            
            function clearMotionFcn(src, ~, motionFcn)
                src.WindowButtonMotionFcn = motionFcn;
                d.layout();
            end
        end
        %------------------------------------------------------------------
        function scrollMenu(d, src, event)
            mw = d.layoutParams.menuWidth;
            p = src.CurrentPoint;
            if p(1) < mw % over menu
                fp = d.Figure.Position;
                mp = d.Menu.Position;
                if event.VerticalScrollCount > 0
                    if mp(2) < 5
                        n = event.VerticalScrollAmount;
                        d.Menu.Position(2) = min(mp(2)+n*30, 5);
                    end
                elseif event.VerticalScrollCount < 0
                    if mp(2)+ mp(4) > fp(4)-10
                        n = event.VerticalScrollAmount;
                        d.Menu.Position(2) = max(mp(2)-n*25, fp(4)-mp(4)-1);
                    end
                end
            end
        end
        % -----------------------------------------------------------------
        function setMousePointer(d, ~, ~)
            [mw, ah, sp] = deal(d.layoutParams.menuWidth, d.layoutParams.distAxisHeight, d.layoutParams.itemSpace);
            p = d.Figure.CurrentPoint;
            if p(1) <= mw + sp
                if p(1) <= mw % inside menu panel
                    set(d.Figure,'Pointer','arrow');
                elseif ~any(structfun(@(x)strcmp(x.Enable, 'on'), d.interactiveModes)) % resize not allowed in interactive mode
                    set(d.Figure,'Pointer','left');
                else
                    set(d.Figure,'Pointer','arrow');
                end
            elseif (p(2)>ah+2*sp) && (p(2)<ah+3*sp) % between axes
                if ~any(structfun(@(x)strcmp(x.Enable, 'on'), d.interactiveModes))
                    set(d.Figure,'Pointer','top');
                else
                    set(d.Figure,'Pointer','arrow');
                end
            elseif ~any(structfun(@(x)strcmp(x.Enable, 'on'), d.interactiveModes))
                set(d.Figure,'Pointer','arrow');
            end
        end
        % -----------------------------------------------------------------
        %% --------- FUNCTIONS PLOTTING -----------------------------------
        function plotPie(d, ax, data, names, cmap)
            [data, pind] = sort(data, 'descend');
            dsum = sum(data);
            ix = find(cumsum(data) >= .95*dsum, 1, 'first');
            [data, pind] = deal(data(ix:-1:1), pind(ix:-1:1));
            h = pie(ax, max([data dsum-sum(data)]/dsum,eps));
            cmap = [cmap(pind,:); .8 .8 .8];
            for i=0:numel(pind)
                set(h(i*2+1),'FaceColor',.85*cmap(i+1,:)+[.15 .15 .15]);
            end
            legend(ax, [names(pind); 'others'], 'Location','EastOutside', 'Interpreter','none');
            set(ax,'PickableParts','all');set(h,'HitTest','off')
        end
        % -----------------------------------------------------------------
        function plotPLT(d, ax, wp, cmap)
            nseg = numel(wp.z);
            if nseg == 1  % Need to trick Matlab
                z = [wp.z; wp.z+1];
                a = [cumsum(wp.alloc,1); zeros(1,numel(wp.alloc))];
                bwidth = 2*numel(wp.z)+3;
                h=barh(ax, z, a,'stacked','BarWidth',bwidth, 'EdgeColor','none');
                args = {'YDir', 'YLim'};
            else
                z = flipud(wp.z);
                if any(diff(z)>=0), z=(nseg:-1:1).'; end
                a = cumsum(flipud(wp.alloc),1);
                if nseg<21
                    args = {'YDir', 'YLim'};
                    h=barh(ax, z, a,'stacked','BarWidth',1, 'EdgeColor','none');
                else
                    h=area(ax, z, a, eps,'EdgeColor','none'); axis(ax,'tight');
                    %                     ax.View(90,-90);
                    args = {'XDir', 'XLim'};
                end
            end
            [zm, zM] = deal(min(z), max(z));
            for i=1:numel(h) % avoid using black bars
                set(h(i),'FaceColor',cmap(i,:))
            end
            set(h,'HitTest','off')
            set(ax,args{1},'reverse', args{2}, [zm zM + sqrt(eps)] + [-.1 .1]*(zM-zm));
        end
        % -----------------------------------------------------------------
        function plotPLT3D(d, ax, inj, reverseY, cmap)
            alloc   = cellfun(@(x) cumsum(x.alloc(end:-1:1,:), 1), inj, 'UniformOutput', false);
            hold(ax,'on');         
            for i = 1:numel(alloc)
                s = abs(sum(alloc{i}));
                mm = max(s);
                pltbarsIx = find(s>0.01*mm);
                if any(pltbarsIx)
                    pltalloc = alloc{i}(:,pltbarsIx);
                    z     = (size(alloc{i},1):-1:1).';
                    h=bar3h(ax, z, pltalloc,'stacked');
                    for j=1:numel(h)
                        set(h(j), 'xdata', get(h(j),'xdata')+i-1, ...
                            'FaceColor', cmap(pltbarsIx(j),:),'EdgeAlpha',.5);
                    end
                    set(h,'HitTest','off')
                end
            end
            box(ax,'on'); grid(ax,'on'); view(ax,[140,20]); axis(ax,'tight');
            set(ax,'Xdir','reverse');
            if reverseY, set(ax,'Ydir','reverse'), end
        end
        % -----------------------------------------------------------------
        function showFPhi(d, ax, modsel, tsel, wsel)
            if isempty(modsel.ix), return, end
            m = getDynamicMeasures(d, modsel, wsel);
            
            singleInj  = numel(wsel.injectorIx)==1;
            singleProd = numel(wsel.producerIx)==1;
            singlePlot = singleInj && singleProd;
            lw = 0.5+singlePlot*1.5; np=0;
            nM = numel(modsel.ix);
            if nM==1
                if singleInj
                    colM = d.prodColors(wsel.producerIx,:);
                elseif singleProd
                    colM = d.injColors(wsel.injectorIx,:);
                else
                    colM = get(gca,'ColorOrder');
                end
                colM = [1 0 0; colM]; lw = 1;
            elseif nM < 8
                colM = get(gca,'ColorOrder');
            else
                colM = colorcube(nM+1);
            end
            hold(ax,'on');
            for n=1:nM
                % Plot F-Phi for the whole regions
                col = colM(n,:);
                if ~singlePlot
                    plot(ax,m.Phi,m.Ft(:,n), ...
                        'LineWidth',2,'Color',col,'DisplayName', modsel.modelNames{modsel.ix(n)});
                end
                if ~m.computePairs, continue; end
                
                % Plot F-Phi for individual wall pairs. For multiple time
                % steps, these lines are plotted using a lighter color. For a
                % single time step, we compute using different colors
                np = size(m.F,2);
                if nM>1
                    if ~singlePlot
                        col = repmat(.5*(col + [.6 .6 .6]), np, 1);
                    end
                elseif np<size(colM,1)
                    col = colM(2:end,:);
                else
                    col = colorcube(np+1);
                end
                for i=1:np
                    plot(ax,m.Phi,m.F(:,i,n),'Color',col(i,:), 'LineWidth',lw,...
                        'DisplayName', m.names{i});
                end
            end
            hold(ax,'off');
            if ~isempty(m.wellName)
                title(ax,['Well: ' m.wellName]);
            end
            if singlePlot
                title(ax,['Well-pair: ' m.wellName ',' m.names{:}]);
                lgn = legend(ax,  modsel.modelNames{modsel.ix});
                set(lgn,'FontSize',8);
            elseif nM<2
                lgn = legend(ax);
                set(lgn,'FontSize',8);
            else
                hax = get(ax,'Children');
                lgn = legend(hax(np+1:np+1:end), ...
                    modsel.modelNames{modsel.ix});
                set(lgn,'FontSize',8);
            end
            ax.XLabel.String = '\phi';
            ax.YLabel.String = 'F';
        end
        % -----------------------------------------------------------------
        function showFracRecovery(d, ax, modsel, tsel, wsel)
            if isempty(modsel.ix), return, end
            m = getDynamicMeasures(d, modsel, wsel);
            
            singleInj  = numel(wsel.injectorIx)==1;
            singleProd = numel(wsel.producerIx)==1;
            singlePlot = singleInj && singleProd;
            lw = 0.5+singlePlot*1.5; np=0;
            nM = numel(modsel.ix);
            if nM==1
                if singleInj
                    colM = d.prodColors(wsel.producerIx,:);
                elseif singleProd
                    colM = d.injColors(wsel.injectorIx,:);
                else
                    colM = get(gca,'ColorOrder');
                end
                colM = [1 0 0; colM]; lw = 1;
            elseif nM < 8
                colM = get(gca,'ColorOrder');
            else
                colM = colorcube(nM+1);
            end
            hold(ax,'on');
            for n=1:nM
                % Plot fractional recovery for the whole regions
                col = colM(n,:);
                if ~singlePlot
                    plot(ax,m.tDt(:,n),1-m.Ft(:,n), ...
                        'LineWidth',2,'Color',col,'DisplayName',  modsel.modelNames{modsel.ix(n)});
                end
                if ~m.computePairs, continue; end
                
                % Plot fractional recovery for individual wall pairs. For multiple time
                % steps, these lines are plotted using a lighter color. For a
                % single time step, we compute using different colors
                hold(ax,'on');
                np = size(m.tD,2);
                if nM>1
                    if ~singlePlot
                        col = repmat(.5*(col + [.6 .6 .6]), np, 1);
                    end
                elseif np<size(colM,1)
                    col = colM(2:end,:);
                else
                    col = colorcube(np+1);
                end
                for i=1:np
                    plot(ax,m.tD(:,i,n),1-m.F(:,i,n),'Color',col(i,:), 'LineWidth',lw, 'DisplayName', m.names{i});
                end
            end
            if ~isempty(m.wellName)
                title(ax,['Well: ' m.wellName]);
            end
            if singlePlot
                title(ax,['Well-pair: ' m.wellName ',' m.names{:}]);
                lgn = legend(ax,  modsel.modelNames{modsel.ix});
                set(lgn,'FontSize',8);
            elseif nM<2
                if np == 0
                    lgn = legend(ax,  modsel.modelNames{modsel.ix});
                else
                    lgn = legend(ax);
                end
                set(lgn,'FontSize',8);
                
            else
                hax = get(ax,'Children');
                lgn = legend(hax(np+1:np+1:end), ...
                    modsel.modelNames{modsel.ix});
                set(lgn,'FontSize',8);
            end
            hold(ax,'off'); axis(ax,'tight');
            set(ax,'XLim',[0 min(5,max(m.tDt(:)))]);
            ax.XLabel.String = 't_D';
            ax.YLabel.String = '1-F';
        end
        % -----------------------------------------------------------------
        function showSweep(d, ax, modsel, tsel, wsel)
            if isempty(modsel.ix), return, end
            m = getDynamicMeasures(d, modsel, wsel);
            
            singleInj  = numel(wsel.injectorIx)==1;
            singleProd = numel(wsel.producerIx)==1;
            singlePlot = singleInj && singleProd;
            lw = 0.5+singlePlot*1.5; np=0;
            nM = numel(modsel.ix);
            if nM==1
                if singleInj
                    colM = d.prodColors(wsel.producerIx,:);
                elseif singleProd
                    colM = d.injColors(wsel.injectorIx,:);
                else
                    colM = get(gca,'ColorOrder');
                end
                colM = [1 0 0; colM]; lw = 1;
            elseif nM < 8
                colM = get(gca,'ColorOrder');
            else
                colM = colorcube(nM+1);
            end
            hold(ax,'on');
            for n=1:nM
                % Plot sweep for the whole regions
                col = colM(n,:);
                if ~singlePlot
                    plot(ax,m.tDt(:,n),m.Ev, ...
                        'LineWidth',2,'Color',col,'DisplayName',  modsel.modelNames{modsel.ix(n)});
                end
                if ~m.computePairs, continue; end
                
                % Plot sweep for individual wall pairs. For multiple time
                % steps, these lines are plotted using a lighter color. For a
                % single time step, we compute using different colors
                hold(ax,'on');
                np = size(m.tD,2);
                if nM>1
                    if ~singlePlot
                        col = repmat(.5*(col + [.6 .6 .6]), np, 1);
                    end
                elseif np<size(colM,1)
                    col = colM(2:end,:);
                else
                    col = colorcube(np+1);
                end
                for i=1:np
                    plot(ax,m.tD(:,i,n),m.Ev,'Color',col(i,:), 'LineWidth',lw, 'DisplayName', m.names{i});
                end
            end
            if ~isempty(m.wellName)
                title(ax,['Well: ' m.wellName]);
            end
            if singlePlot
                title(ax,['Well-pair: ' m.wellName ',' m.names{:}]);
                lgn = legend(ax,  modsel.modelNames{modsel.ix});
                set(lgn,'FontSize',8);
            elseif nM<2
                if np == 0
                    lgn = legend(ax,  modsel.modelNames{modsel.ix});
                else
                    lgn = legend(ax);
                end
                set(lgn,'FontSize',8);
                
            else
                hax = get(ax,'Children');
                lgn = legend(hax(np+1:np+1:end), ...
                    modsel.modelNames{modsel.ix});
                set(lgn,'FontSize',8);
            end
            hold(ax,'off'); axis(ax,'tight');
            set(ax,'XLim',[0 min(5,max(m.tDt(:)))]);
            ax.XLabel.String = 't_D';
            ax.YLabel.String = 'E_v';
        end
        % -----------------------------------------------------------------
        function showLorenz(d, ax, modsel, tsel, wsel)
            if isempty(modsel.ix), return, end
            m = getDynamicMeasures(d, modsel, wsel);
            
            singlePlot  = (numel(wsel.injectorIx)==1) && ...
                (numel(wsel.producerIx)==1);
            axes(ax)
            nM = numel(modsel.ix);
            if ~m.computePairs
                h=bar(ax, m.LCt,'LineWidth',1, ...
                    'FaceColor',[.9 .9 .9],'EdgeColor',[0 .5 .8]);
                hold(ax,'on');
                text(h.XData, h.YData, cellstr(num2str(m.LCt(:),'%.4f')), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment','bottom', 'FontSize', 8);
                set(ax,'XTickLabelRotation',45);
                ax.YLabel.String = 'L_c';
                
            elseif singlePlot
                h=bar(ax, m.LC,'LineWidth',1, ...
                    'FaceColor',[.9 .9 .9],'EdgeColor',[0 .5 .8]);
                hold(ax,'on');
                text(h.XData, h.YData, cellstr(num2str(m.LC(:),'%.4f')), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment','bottom', 'FontSize', 8);
                set(ax,'XTickLabel', modsel.modelNames(modsel.ix));
                set(ax,'XTickLabelRotation',45,'FontSize',8);
                ax.YLabel.String = 'L_c';
            elseif nM<4
                n = size(m.LC,1);
                h=bar(ax, m.LC,'LineWidth',1, ...
                    'FaceColor',[.9 .9 .9],'EdgeColor',[0 .5 .8]);
                grid(ax,'on'); hold(ax,'on');
                plot(ax,[.65 n+.35],[m.LCt m.LCt]','Color','r', 'LineWidth',2);
                if nM==1
                    text(h.XData, h.YData, cellstr(num2str(m.LC(:),'%.4f')), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment','bottom', 'FontSize', 8);
                else
                    yb = cat(1, h.YData);
                    xb = bsxfun(@plus, h(1).XData, [h.XOffset]');
                    text(xb(:),yb(:)-.05, cellstr(num2str(yb(:),'%.4f')), ...
                        'rotation', 90, 'horiz', 'right', 'FontSize',6);
                end
                set(ax,'XTickLabel',m.names,'FontSize',8);
                if n>8, set(ax,'XAxisLocation','top','XTickLabelRotation',45), end
                ax.YLabel.String = 'L_c';
            else
                n = size(m.LC,1);
                bar3(ax,m.LC); view(-80,20)
                hold(ax,'on')
                plot3(repmat(1:nM,2,1), repmat([.65 n+.35],nM,1)', ...
                    [m.LCt m.LCt]','Color','r','LineWidth',2);
                set(ax,'YTickLabel',m.names,'FontSize',8)
                if numel(m.names)>8, set(ax,'YTickLabelRotation',45), end
                ax.ZLabel.String = 'L_c';
            end
            if ~isempty(m.wellName), title(['Well: ' m.wellName]); end
            hold(ax,'off'); axis(ax,'on', 'normal');
            
        end
        % -----------------------------------------------------------------
        function PlotTOFArrival(d, ax, modsel, tsel, wellIx, wellType, extendTime)
            
            m = modsel.ix;
            PORVIx = arrayfun(@(x) strcmp(x.name,'PORV'),d.Data{m}.static);
            pv = d.Data{m}.static(PORVIx).values;
            D = d.Data{m}.diagnostics.D;
            state = d.Data{m}.states;
            
            [data,tof] = getPhaseTOFDistributions(state, pv, wellIx, wellType, [], [], D);
            h = area(ax, tof, data);

            phcol = {[.4 .4 1],[.3 0 0], [.5 1 .5]};
            for ph = 1:numel(h)
                set(h(ph), 'FaceColor', phcol{ph});
            end
            
            %axis([min(tof_sub), max(tof_sub), 0, 1])
            hold(ax,'off'); axis(ax,'tight');
            ax.XLabel.String = 'TOF distance in years';
            
            phnames = {d.Data{1}.dynamic.name};
            legend(ax, phnames(2:end), 'Location', 'Northwest');
            set(ax, 'XLim', [0, min(ax.XLim(2),extendTime)]);
        end
        
        % -----------------------------------------------------------------
        function PlotPhaseTOFArrivals(d, ax, modsel, tsel, wsel, wtype, phase, extendTime)
            m      = modsel.ix;
            PORVIx = arrayfun(@(x) strcmp(x.name,'PORV'),d.Data{m}.static);
            pv     = d.Data{m}.static(PORVIx).values;
            D      = d.Data{m}.diagnostics.D;
            state  = d.Data{m}.states;
            switch wtype
                case 'producer'
                    wellIx1   = wsel.producerIx;
                    wellIx2   = wsel.injectorIx;
                    com       = wsel.communicationMatrix;
                    colors    = d.injColors;
                    wnames   = d.WellPlot.injectorNames;
                case 'injector'
                    wellIx1   = wsel.injectorIx;
                    wellIx2   = wsel.producerIx;
                    com       = wsel.communicationMatrix';
                    colors    = d.prodColors;
                    wnames   = d.WellPlot.producerNames;
                otherwise
                    error('Incorrect well type');
            end
                      
            % If communicating wells are not selected: autodetect well-pairs
            if isempty(wellIx2)
                tot = sum(com(:));
                n   = max(size(com));
                com = com > (wsel.communicationLimit/100)*tot/n;
                ix = logical(sum(com(:,wellIx1), 2));
                wellIx2 = find(ix).';
            end
            
            % Compute phase distributions of individual phases
            [sphase,tof] = ...
                getPhaseTOFDistributions(state, pv, wellIx1, wtype, phase, wellIx2, D);
            
            % Add 0 to start to make plot look better
            h = area(ax, [0; tof], [zeros(1,size(sphase,2)); sphase]);
            for ph = 1:numel(h)
                set(h(ph), 'FaceColor', colors(wellIx2(ph),:));
            end
            hold(ax,'off'); axis(ax,'tight');
            ax.XLabel.String = 'TOF distance in years';
            ax.YLabel.String = 'm^3';
            set(ax, 'XLim', [0, min(ax.XLim(2),extendTime)]);
            ymax = 0;
            for i  = 1:numel(h)
                ymax = ymax + interp1(h(i).XData, h(i).YData, ax.XLim(2));
            end
            if ymax > 0
                set(ax, 'YLim', [0, ymax+ymax.*0.01]);
            end
            legend(ax, wnames(wellIx2), 'Location', 'Northwest');
        end
        
        % -----------------------------------------------------------------
        function infoStr = getModelInfoString(d,s3)
            m = getSelectedModelIndex(d,s3);           
            plotProp3D = s3.psel.propPopup.String{s3.psel.propIx};
            plotType3D = s3.psel.typePopup.String{s3.psel.typeIx};
            if numel(m) > 1
                plotStat3D = s3.psel.statPopup.String{s3.psel.statIx};
            else
                plotStat3D = "n/a";
            end
            
            infoStr1 = sprintf("Displayed Models:\n");
            infoStr2 = sprintf('    %s\n',s3.modsel.modelNames{m});
            infoStr3 = sprintf("\nDisplayed Property: \n    %s (%s) \nStatistic: \n    %s", ...
                plotProp3D, plotType3D, plotStat3D);
            infoStr = strcat(infoStr1,infoStr2,infoStr3);
            
        end
        
        
        %% --------- FUNCTIONS GENERAL ------------------------------------
        function resetSelectors(d, s, selector, field)
            flds = reshape(fieldnames(s),1, []);
            for f = flds(~strcmp(flds, selector))
                if isprop(s.(f{1}), field)
                    s.(f{1}).(field) = 1;
                end
            end
        end
        % -----------------------------------------------------------------
        function setColormap3D(d, str)
            if ~strcmp(str, d.colormap3D)
                if strcmp(str, 'injectors')
                    ninj = d.WellPlot.nInj;
                    cmap = d.injColors(1:ninj,:);
                    colormap(d.Axes3D, cmap);
                elseif strcmp(str, 'producers')
                    nprod = d.WellPlot.nProd;
                    cmap = d.prodColors(1:nprod,:);
                    colormap(d.Axes3D, cmap);
                else
                    colormap(d.Axes3D, str);
                end
                d.colormap3D = str;
            end
        end
        
        % -----------------------------------------------------------------
        function updateColorHist(d)
            vals = d.Patch.colorData(d.Patch.cells);
            lims = d.Axes3D.CLim;
            [counts,bins] = hist(vals, linspace(lims(1), lims(2),26));
            visState = 'on';
            c = d.colorHAx.Children;
            if numel(c) == 1 && isprop(c, 'Visible')
                visState = d.colorHAx.Children.Visible;
            end
            barh(d.colorHAx, bins, counts,'hist');
            axis(d.colorHAx,'tight', 'off');
            h = findobj(d.colorHAx,'Type','patch');
            set(h,'FaceColor',[.8 .9 1],'EdgeColor',[.2 .2 .3], 'Visible', visState);
        end
        % -----------------------------------------------------------------
        function updateColorBar(d, s3, lims, islog)
            if nargin < 4
                islog = false;
            end
            cb = d.colorBar;
            % let first matlab choose ticks/labels
            if ~(s3.psel.typeIx == 3 && any(s3.psel.propIx == [7,8]))
                cb.TicksMode = 'auto';
                cb.TickLabelsMode = 'auto';
                cb.Limits = lims;
                if islog
                    p = cb.Ticks;
                    cb.TickLabels = cellfun(@(x)sprintf('%2.1e',10^x), num2cell(p), ...
                        'UniformOutput', false);
                end
                set(d.colorHAx.Children, 'Visible', 'on');
                d.setColormap3D('default');
                %cb.Position(1) = d.layoutParams.menuWidth + 50;
                %d.colorHAx.Position(1) = d.colorBar.Position(1)+d.colorBar.Position(3)+10;
            else
                if  s3.psel.propIx == 7 % sweep regions
                    ninj = d.WellPlot.nInj;
                    d.setColormap3D('injectors');
                    cb.Ticks = (.5:ninj)/ninj;
                    cb.TickLabels = s3.wsel.injSelector.String;
                    cb.Limits = [0 1];
                elseif s3.psel.propIx == 8 % drainage regions
                    nprod = d.WellPlot.nProd;
                    d.setColormap3D('producers');
                    cb.Ticks = (.5:nprod)/nprod;
                    cb.TickLabels = s3.wsel.prodSelector.String;
                    cb.Limits = [0 1];
                end
                set(d.colorHAx.Children, 'Visible', 'off');
                cb.TickLabelInterpreter = 'none';
                %cb.Position(1) = d.layoutParams.menuWidth + 75;
                % switch off log if selected
                if s3.psel.logSwitch == true
                    s3.psel.logSwitch = false;
                    s3.psel.logSwitchCallback(s3.psel.logSwitchBox, []);
                end
            end
            d.colorBar = cb;
        end
        % -----------------------------------------------------------------
        
        function m = getSelectedModelIndex(d,s3)
            
            if isempty(s3.modsel.ix)
                m = 1;
                s3.modInfo.infoText = (s3.modsel.modelNames{m});
            else
                m = s3.modsel.ix;
            end
        end
        
        % -----------------------------------------------------------------    
        function wellCommunication = getGlobalCommunicationMatrix(d)
            % Finds max connection between each well pair for all
            % models.
            nM = numel(d.Data);
            allCom = zeros([size(d.Data{1}.diagnostics.wellCommunication) nM]);
            for i = 1:nM
                allCom(:,:,i) = d.Data{i}.diagnostics.wellCommunication;
            end
            
            wellCommunication = max(allCom,[],3);
        end
        
        % -----------------------------------------------------------------
        function wellCommunication = getSelectedModelCommunicationMatrix(d,s3)
            % Finds max connection between each well pair for all selected
            % models.
            if isempty(s3.modsel.ix)
                return
            else
                m = s3.modsel.ix;
                nM = numel(m);
                allCom = zeros([size(d.Data{1}.diagnostics.wellCommunication) nM]);
                for i = 1:nM
                    allCom(:,:,i) = d.Data{m(i)}.diagnostics.wellCommunication;
                end
                wellCommunication = max(allCom,[],3);
             end
        end
        % -----------------------------------------------------------------
        
        function d = switchTransparancy(d, st)
            if strcmp(st, 'on')
                alp = [.2, .5];
            else
                alp = [1, 1];
            end
            d.outlineGrid.EdgeAlpha = alp(1);
            d.Patch.patchMain.EdgeAlpha = alp(2);
            ii = find(strcmp(d.Patch.patchOpt, 'EdgeAlpha'));
            if numel(ii) == 1
                d.Patch.patchOpt{ii+1} = alp(2);
            end
        end
        
        
    end
end % --- END POSTPROCESSDIAGNOSICS ---------------------------------------

%% -- UTILITY FUNCTIONS ---------------------------------------------------
function [models,wells,states0] = extractInput(arg1, arg2)
if iscell(arg1) && iscell(arg2)
    models  = arg1;
    wells   = arg2;
    states0 = cell(size(arg1));
elseif isa(arg1,'ModelEnsemble') && isnumeric(arg2)
    [models, wells, states0] = deal(cell(1,numel(arg2)));
    fprintf(1,'Starting to extract ensemble members\n');
    for k = 1:numel(arg2)
        fprintf(1,'  member no %d ... ', arg2(k));
        tmp       = arg1.setupFn(arg2(k));
        models{k} = tmp.model;
        wells{k}  = tmp.W;
        if isfield(tmp, 'state0')
            states0{k} = tmp.state0;
        end
        fprintf(1,'done\n');
    end
    fprintf(1,'Extracted %d ensemble members\n', numel(arg2));
else
    error('Incorrect input specification of model ensemble');
end
end
% -----------------------------------------------------------------
function s = emptyDiagnostics(d)
s = struct('name', '', 'values', []);
str = d.Props.diagnostics.name;
s = repmat(s, 1, numel(str));
for k = 1:numel(str)
    s(k).name = str{k};
end
end
% -----------------------------------------------------------------
function [lims, vals, flag] = makeSafeForLog(lims, vals, oom)
if lims(2) <= 0
    flag = false;
else
    flag = true;
    if lims(1) <= 0 || lims(1) < lims(2)*10^(-oom)
        lims(1) = lims(2)*10^(-oom);
        vals = max(vals, lims(1));
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
function varargout = addPatchContextMenu(varargin)
nax = nargin;
varargout = cell(1,nax);
for k = 1:nax
    input = varargin{k};
    inputIsPatch = true;
    p = [];
    if isa(input, 'matlab.graphics.primitive.Patch')
        p = input;
    elseif isa(input, 'CellDataPatch')
        p = input.patchMain;
        inputIsPatch = false;
    end
    
    varargout{k} = input;
    
    if ~isempty(p)
        m = uicontextmenu('Parent', p.Parent.Parent);
        uimenu(m, 'Label', 'Export to new figure', 'Callback', @(src, event)copyToFig(src, event, p.Parent))
        if inputIsPatch
            varargout{k}.UIContextMenu = m;
        else
            varargout{k}.patchMain.UIContextMenu = m;
        end
    else
        warning('Unexpteted input: %s' ,class(input));
    end
end
end
% -----------------------------------------------------------------
function copyToFig(src, event, ax) 
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
% -----------------------------------------------------------------
% Toolbar customization
function d = addExtraTools(d)
tb  = findall(d.Figure,'Type','uitoolbar');
gh = uitoggletool(tb,'CData', getIcon('grid_outline.png'), 'Separator', 'on', 'State', 'on', ...
    'TooltipString', 'Show/hide outline grid');
gh.OnCallback = @(src, event)outlineOnOff(src, event, d, 'on');
gh.OffCallback = @(src, event)outlineOnOff(src, event, d, 'off');

lh = uitoggletool(tb,'CData',getIcon('light.png'), 'Separator', 'off', 'State', 'off', ...
    'TooltipString', 'Enable/disable light');
lh.OnCallback = @(src, event)lightOnOff(src, event, d, 'on');
lh.OffCallback = @(src, event)lightOnOff(src, event, d, 'off');

lh = uitoggletool(tb,'CData',getIcon('alpha.png'), 'Separator', 'on', 'State', 'on', ...
                   'TooltipString', 'Enable/disable transparancy (disable for faster rendering)');
lh.OnCallback = @(src, event)alphaOnOff(src, event, d, 'on');
lh.OffCallback = @(src, event)alphaOnOff(src, event, d, 'off');

aspFac = 1.2;
zp = uipushtool(tb,'CData',getIcon('z_plus.png'), 'Separator', 'on', ...
    'TooltipString','Increase vertical aspect ratio');
zp.ClickedCallback = @(src, event)changeVerticalAspectRatio(src, event, d, 1./aspFac);

zm = uipushtool(tb,'CData',getIcon('z_minus.png'), 'Separator', 'on', ...
    'TooltipString','Decrease vertical aspect ratio');
zm.ClickedCallback = @(src, event)changeVerticalAspectRatio(src, event, d, aspFac);

end
function outlineOnOff(~, ~, d, st)
d.outlineGrid.Visible = st;
end

function alphaOnOff(~,~,d,st)
d.switchTransparancy(st);
end

function lightOnOff(~, ~, d, st)
if isvalid(d.camlight)
    d.camlight.Visible = st;
else
    d.camlight = camlight;
end
end
function changeVerticalAspectRatio(~, ~, d, fac)
d.Axes3D.DataAspectRatio(3) = fac*d.Axes3D.DataAspectRatio(3);
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
