classdef TrajectoryGUI < handle
    properties
        model
        state0
        wellCases
        Data
        % graphics
        Figure
        Axes                = gobjects; % (1)left axis, (2) right axis, (3) bottom left  
        WellPlot1
        WellPlot2
        Menu
        Patch3D
        subGrid             % grid for subregion            
        ccol                % colordata for patch
        gradientVectors     = gobjects;
        % special graphics
        xyRegion            % interactive xy-subregion
        lineOnSlice         % interactive trajectory with slice    
        trajInMain          % 'copy' of trajectory in Axes(1)
        trajInMainPnt       % interactive point in Axes(1)
        % selectors
        caseSelector
        wellSelector
        sliceSelector
        propertySelector
        proxySelector
        %
        interactiveModes    % struct with rotate/zoom/pan - options        
        layoutParams   = struct('menuWidth',        300, ...
                                'itemSpace',         40, ...
                                'includeToolTips', true);
        wellAuxiliaries     % aux well-view/plotting variables
        % proxy-related
        obj
        objValues
    end
    properties (Dependent)
        caseNo
        wellNo
    end
    
    methods
        function d = TrajectoryGUI(model, W, varargin)
            require wellpaths
            opt = struct('style', 'default', 'state0', [], 'objective', []);
            opt = merge_options(opt, varargin{:});
            
            % ------ Create figure window ---------------------------------
            screensize = get( groot, 'Screensize' );
            wsize = .75*screensize(3:4);
            wsize = [screensize(3:4)-wsize*1.1, wsize];
            d.Figure = limitedToolbarFigure('Position', wsize);
            
            d.model  = model;
            d.state0 = opt.state0;
            d.wellCases{1} = addTrajectories(W, d.model.G, 10);
            d.wellAuxiliaries = repmat(struct('nPoints', [], 'rPos', []), ...
                                       [numel(W), 1]);
            % add fields for faster comp
            d.model.G = addBoundingBoxFields(d.model.G);
            % check objective
            hasObj = false;
            if ~isempty(opt.objective)
                if isa(opt.objective, 'DiagnosticsObjective')
                    d.obj = opt.objective;
                    hasObj = true;
                else
                    warning('Objective of unsupported type, disregarded');
                end
            end
            % get cell-data
            d.Data = getCellData(model, d.state0);
            props  = struct('static', struct('name', {{d.Data.name}}, ...
                                             'limits', {{d.Data.limits}}));
            %--------------------------------------------------------------
            itemOpts = {'Parent', d.Figure, 'Visible','off', 'style', opt.style};
            
            csel = WellCaseSelector(itemOpts{:});
            wsel = WellEditSelectorNew(d.wellCases{1}, 'wellNames', ...
                {d.wellCases{1}.name}, itemOpts{:});
            ssel = TrajectorySliceSelector(itemOpts{:});
            psel = PropertyDisplaySelector('props', props, itemOpts{:}, ...
                        'includeLogSwitch', true, ...
                        'includeFilter',    true);
            xsel = ProxySelector('Parent', d.Figure);      
            
            d.Menu = UIMenu('Title', 'Menu', 'Parent', d.Figure, ...
                itemOpts{:}, 'items', {csel, wsel, ssel, psel, xsel});
            
            tt = [251 180 76; 58 152 216; 42 187 155; 252 121 122; 155, 89, 182]./255;
            tt = [tt; tatarizeMap];
            for k =1:numel(d.Menu.items)
                d.Menu.items{k}.BackgroundColor = tt(k,:);
                try
                    d.Menu.items{k}.titleColor = 'w';
                catch
                    d.Menu.items{k}.ForegroundColor = 'w';
                end
            end
            % setup axes
            d.Axes(1) = axes('Parent', d.Figure, 'Units', 'pixels', 'ZDir', 'reverse');
            d.Axes(2) = axes('Parent', d.Figure, 'Units', 'pixels', 'ZDir', 'reverse');
            d.ccol  = model.rock.perm(:,1);
            d.Patch3D = CellDataPatch(model.G, d.ccol, ...
                'Parent', d.Axes(2), 'EdgeColor', [.3 .3 .3], ...
                'FaceAlpha', .2, 'EdgeAlpha', .3, ...
                'BackFaceLighting', 'lit', ...
                'Hittest', 'off');
            d.updateWellPlots();
            % setup interactive modes
            d.interactiveModes = setCustomInteractiveModes(d.Figure);
            d.interactiveModes.rotate3d.setAllowAxesRotate(d.Axes(2), false);
            % Callbacks
            wsel.changedCallback = @d.enableSelectorsCallback;
            wsel.Callback        = @d.wellSelectorCallback;
            csel.Callback        = @d.caseSelectorCallback;
            ssel.Callback        = @d.sliceControlCallback;
            psel.Callback        = @d.propertySelectorCallback;
            xsel.Callback        = @d.proxySelectorCallback;
            
            d.Figure.SizeChangedFcn        = @d.layout;
            d.layout();
            
            wsel.typePopup.Enable = 'off';
            d.wellSelector  = wsel;
            d.caseSelector  = csel;
            d.sliceSelector = ssel;
            d.propertySelector = psel;
            d.proxySelector = xsel;
            
            if ~hasObj
                d.proxySelector.Visible = 'off';
            else
                d.Axes(3) = axes('Parent', d.Figure, 'Units', 'pixels');
            end
            
            % run callback for first selection
            d.wellSelectorCallback();
            % set various graphics
            axis(d.Axes(1), 'vis3d');
            view(d.Axes(1), 3);
            daspect(d.Axes(1), [1 1 .4]);
            zoom(d.Axes(1), .7)
            zoom(d.Axes(1), 'reset')
            set(d.Axes(1), {'XLimMode', 'YLimMode', 'ZLimMode'}, {'manual','manual','manual'});
            zl = d.Axes(1).ZLim;
            d.Axes(1).ZLim = zl + [-.15 .15]*diff(zl, [], 2);
            % select PERMX for display
            d.propertySelector.propIx = 2;
            d.propertySelector.Callback();

        end
        
        function set.wellNo(d, val)
            d.wellSelector.wellNo = val;
        end
        function val = get.wellNo(d)
            val = d.wellSelector.wellNo;
        end
        
        function set.caseNo(d, val)
            d.caseSelector.ix = val;
        end
        function val = get.caseNo(d)
            val = d.caseSelector.ix;
        end
        
       %-------------------------------------------------------------------
       % Don't allow case-selection until well-changes are saved/reset
        function enableSelectorsCallback(d, src, event)
            if strcmp(d.wellSelector.enableSaveReset, 'on')
                d.caseSelector.Enable = 'off';
            else
               d.caseSelector.Enable = 'all';
            end
        end
        
        %-----------------------------------------------------------------
        function caseSelectorCallback(d, src, event)
            switch src.Tag
                case 'select'
                    d.wellSelector.updateWells(d.wellCases{d.caseNo(1)});
                    d.updateWellPlots();
                    d.WellPlot1.visibleCases = d.caseNo;
                    d.WellPlot2.visibleCases = d.caseNo;
                case 'new'
                    % copy current case
                    ncase = numel(d.wellCases);
                    d.wellCases = d.wellCases([1:ncase, d.caseNo]);
                    d.wellSelector.updateWells(d.wellCases{d.caseNo});
                    d.updateWellPlots();
                    d.WellPlot1.visibleCases = d.caseNo;
                    d.WellPlot2.visibleCases = d.caseNo;
                    %d.wellAuxiliaries = d.wellAuxiliaries([1:ncase, d.caseNo]);
                case 'delete'
                    d.wellCases(d.caseNo) = [];
                    d.objValues(intersect(d.caseNo, 1:numel(d.objValues))) = [];
                    d.updateWellPlots();
                    d.WellPlot1.visibleCases = 1;
                    d.WellPlot2.visibleCases = 1;
                case 'launch'
                    casenm = d.caseSelector.listbox.String(d.caseNo);
                    DiagnosticsViewer(repmat({d.model}, [1, numel(casenm)]), ...
                        d.wellCases(d.caseNo), 'modelNames', casenm, 'includeAverage', false, ...
                        'state0', repmat({d.state0}, [1, numel(casenm)]));

            end
            %delete(d.gradientVectors);
        end
        
        %-----------------------------------------------------------------
        function wellSelectorCallback(d, src, event)
            if nargin == 1
                % accept empty input
                [src, event] = deal(d.wellSelector.wellPopup, []);
            end
            switch src.Tag
                case 'name'
                    d.updateXYRegion(src, event);
                    %delete(d.gradientVectors);
                case 'save'
                    isNew = d.wellNo > numel(d.wellCases{d.caseNo});
                    % update trajectory
                    tmp  = d.wellSelector.W(d.wellNo);
                    traj = d.lineOnSlice.trajectory;
                    tmp  = updateWellTrajectory(d.model, tmp, traj);
                    d.wellCases{d.caseNo}(d.wellNo) = tmp;
                    d.wellSelector.W(d.wellNo) = tmp;
                    if isNew % add well to other cases and set shut
                        for kc = 1:numel(d.wellCases)
                            if kc ~= d.caseNo % set to shut
                                d.wellCases{kc}(d.wellNo) = d.wellSelector.W(d.wellNo);
                                d.wellCases{kc}(d.wellNo).status = false;
                            end
                        end
                        d.wellAuxiliaries(d.wellNo) = struct('nPoints', d.sliceSelector.nPnt, ...
                                                             'rPos',    d.xyRegion.Position);   
                    end
                    d.updateWellPlots();
                    d.WellPlot1.visibleCases = d.caseNo;
                    d.WellPlot2.visibleCases = d.caseNo;
                    d.wellSelector.enableSaveReset = 'off';
                    % remove previously computed objectives
                    d.objValues{d.caseNo} = [];
            end
        end
        
        %-----------------------------------------------------------------
        function updateXYRegion(d, src, event)
            if strcmp(src.Style, 'popupmenu')
                [cno, wno] = deal(d.caseNo, d.wellNo);
                if isempty(d.wellAuxiliaries(wno).rPos) 
                    c = d.model.G.cells.centroids(d.wellCases{cno(1)}(wno).cells,:);
                    dd = 100;
                    [minp, maxp] = deal(min(c(:,1:2))-dd, max(c(:,1:2))+dd);
                    d.wellAuxiliaries(wno).rPos = [minp, maxp-minp];
                end
                pos = d.wellAuxiliaries(wno).rPos;
                if ~isempty(d.xyRegion)
                    d.xyRegion.Position = pos;
                else
                    d.xyRegion = InteractiveRectangle('Position', pos, 'Parent', d.Axes(2));
                    d.xyRegion.Callback = @d.updateTrajectory;
                end
            end
            d.updateTrajectory(src,event);
            d.updateTrajInMainPnt(src,event);
            d.wellSelector.enableSaveReset = 'off';
            %
        end
        
        %-----------------------------------------------------------------
        function updateTrajectory(d, src, event)
            pos = d.xyRegion.Position;
            d.wellAuxiliaries(d.wellNo).rPos = pos;
            [xlim, ylim] = deal([pos(1), pos(1)+pos(3)], [pos(2), pos(2)+pos(4)]);
            cc = d.model.G.cells.centroids;
            cix = (cc(:,1) >= xlim(1) & cc(:,1) <= xlim(2)) &...
                (cc(:,2) >= ylim(1) & cc(:,2) <= ylim(2));
            d.subGrid = extractSubgrid(d.model.G, cix);
            set(d.Axes(1), {'XLim', 'YLim'}, {xlim, ylim});
            zoom(d.Figure, 'reset');
            % end
            if ~isa(src, 'matlab.ui.Figure') || isempty(d.lineOnSlice)% redraw trajectory
                if ~isempty(d.lineOnSlice) && isvalid(d.lineOnSlice)
                    d.switchInteractionOff();
                    delete(d.lineOnSlice);
                end
                if numel(d.caseNo) == 1
                    coord = d.wellCases{d.caseNo}(d.wellNo).trajectory;
                    d.lineOnSlice = InteractiveLineOnSurf(d.subGrid, 'XData', coord(:,1), ...
                        'YData', coord(:,2), 'ZData', coord(:,3), 'Parent', d.Axes(1), ...
                        'colorData', d.ccol(cix));
                    if isempty(d.wellAuxiliaries(d.wellNo).nPoints)
                        d.wellAuxiliaries(d.wellNo).nPoints = ...
                            d.sliceSelector.nPnt;
                    else
                        d.sliceSelector.nPnt = d.wellAuxiliaries(d.wellNo).nPoints;
                    end
                    d.sliceControlCallback();
                    d.lineOnSlice.Callback = @d.updateTrajInMainPnt;
                end
            elseif isvalid(d.lineOnSlice) % just update grid and slice
                d.lineOnSlice.G = d.subGrid;
                d.lineOnSlice.colorData = d.ccol(cix);
                d.lineOnSlice.updateSlice();
            end
            
        end
        %-----------------------------------------------------------------
        function sliceControlCallback(d, src, event)
            ssel = d.sliceSelector;
            if nargin < 2 %update all
                d.lineOnSlice.setAngle(ssel.angle);
                d.lineOnSlice.redistributePoints(ssel.nPnt);
                d.lineOnSlice.slice.FaceAlpha = ssel.alpha;
                d.lineOnSlice.slice.EdgeAlpha = ssel.alpha;
            else
                if ~isempty(d.lineOnSlice)
                    switch src.Tag
                        case 'angle'
                            d.lineOnSlice.setAngle(ssel.angle);
                        case 'nPnt'
                            d.lineOnSlice.redistributePoints(ssel.nPnt);
                            d.wellAuxiliaries(d.wellNo).nPoints = ssel.nPnt;
                        case 'alpha'
                            d.lineOnSlice.slice.FaceAlpha = ssel.alpha;
                            d.lineOnSlice.slice.EdgeAlpha = ssel.alpha;
                    end
                end
            end
        end
        %-----------------------------------------------------------------
        function propertySelectorCallback(d, src, event)
            pix = d.propertySelector.propIx;
            [v, lims] = deal(d.Data(pix).values, d.Data(pix).limits);
            subix = d.propertySelector.minValue <= v & ...
                    d.propertySelector.maxValue >= v;
            d.Patch3D.cells = subix;
            if d.propertySelector.logSwitch
                [lims, v, flag] = makeSafeForLog(lims, v, 5);
                if ~flag % values not good for log-plot, reset switch
                    d.propertySelector.logSwitchBox.Value = 0;
                else
                    [lims, v] = deal(log10(lims), log10(v));
                end
            end
            d.ccol = v;
            d.Patch3D.colorData = v;
            % For the 2D slice do thresholding only on faces (not edges)
            tmp = v(d.subGrid.cells.global);
            ix  = subix(d.subGrid.cells.global);
            tmp(~ix) = nan;
            d.lineOnSlice.colorData = tmp;
            d.lineOnSlice.updateSlice();
            set([d.Axes(1), d.Axes(2)], 'CLim', lims);
        end
        
        %-----------------------------------------------------------------
        
        function proxySelectorCallback(d, src, event)
            if nargin == 1
                % interpret this as case has been updated -> delete
                % computed proxy-info
                d.objValues{d.caseNo} = [];
                d.proxySelector.slider.Enable = 'off';
                delete(d.gradientVectors);
            else
                switch src.Tag
                    case 'objective'
                        % compute any uncomputed objectives
                        cix = 1:numel(d.caseSelector.names);
                        for k = cix
                            if numel(d.objValues) < k || isempty(d.objValues{k})
                                % reset wellmodels in case number of wells have changed
                                d.obj.model.parentModel.parentModel.FacilityModel.WellModels = {};
                                d.objValues{k} = d.obj.compute(d.wellCases{k}, 'computeGradient', false);
                            end
                        end
                        % plot resulting values
                        plotObjectiveValues(d.Axes(3), d.objValues, d.caseSelector.listbox.String);
                    case {'gradient', 'control'}
                        if numel(d.caseNo) > 1
                            fprintf('Please select single case for gradient computations');
                            return
                        end
                        % compute if needed
                        try 
                            hasGrad = ~isempty(d.objValues{d.caseNo}.gradient.position{d.wellNo});
                        catch
                            hasGrad = false;
                        end
                        if ~hasGrad
                            W = d.wellCases{d.caseNo};
                            W = rmfield(W, 'trajectory'); 
                            W(d.wellNo).posControl = setupPositionControl(d.model.G, W(d.wellNo), d.lineOnSlice);
                            d.obj.model.parentModel.parentModel.FacilityModel.WellModels = {};
                            if ~(numel(d.objValues) < d.caseNo || isempty(d.objValues{d.caseNo}))
                                tmp = d.obj.compute(W, 'computeGradient', true, ...
                                    'state', d.objValues{d.caseNo}.state, ...
                                    'D',     d.objValues{d.caseNo}.D);
                                g = tmp.gradient.position{d.wellNo};
                                d.objValues{d.caseNo}.gradient.position{d.wellNo} = g;
                                d.objValues{d.caseNo}.gradient.well = tmp.gradient.well;
                            else
                                d.objValues{d.caseNo} = d.obj.compute(W, 'computeGradient', true);
                                g = d.objValues{d.caseNo}.gradient.position{d.wellNo};
                            end
                            % scale grad according to axis limits
                            d.objValues{d.caseNo}.gradPlotScale{d.wellNo} = ...
                                getGradPlotScaling(g, get(d.Axes(1), {'XLim', 'YLim', 'ZLim'}));
                        end
                        if strcmp(src.Tag, 'gradient')
                            d.showGradientVectors();
                            d.proxySelector.slider.Enable = 'on';
                        else
                            plotControlGradients(d.Axes(3), d.wellCases{d.caseNo}, d.objValues{d.caseNo}.gradient.well)
                        end
                    case 'slider'
                        d.moveTrajAlongGradient(d.proxySelector.slider.Value);                        
                        
                end
            end
        end
        %-----------------------------------------------------------------
        function updateTrajInMainPnt(d, src, event)
            if ~isempty(d.lineOnSlice) && isvalid(d.lineOnSlice)
                l = d.lineOnSlice.trajectory;
                [x,y] = deal(l.XData(1), l.YData(1));
                if isempty(d.trajInMainPnt) && ~isempty(x)
                    d.trajInMainPnt = InteractivePoint('Parent', d.Axes(2), 'XData', x, 'YData', y, ...
                        'Callback', @d.moveWellXY, ...
                        'moveCallback', @d.updateTrajInMain);
                else
                    d.trajInMainPnt.Position = [x,y];
                end
                d.updateTrajInMain(src, event);
                d.wellSelector.enableSaveReset = 'on';
            end
        end
        %-----------------------------------------------------------------
        function updateWellPlots(d)
            d.WellPlot1 = WellPlotHandle(d.model.G, d.wellCases, 'Parent', d.Axes(1));
            d.WellPlot2 = WellPlotHandle(d.model.G, d.wellCases, 'Parent', d.Axes(2));
            d.WellPlot1.closedColor = 'none';
            d.WellPlot2.closedColor = 'none';
            d.WellPlot1 = setSpecialWellPlotProps(d.WellPlot1);
        end
        %-----------------------------------------------------------------
        function updateTrajInMain(d, src, event)
            if ~isempty(d.lineOnSlice)
                l = d.lineOnSlice.trajectory;
                [x,y,z] = deal(l.XData, l.YData, l.ZData);
                p = d.trajInMainPnt.Position;
                x = x + (p(1)-x(1));
                y = y + (p(2)-y(1));
                if isempty(d.trajInMain) && ~isempty(x)
                    d.trajInMain = line('Parent', d.Axes(2), 'XData', x, 'YData', y, 'ZData', z, ...
                        'Color', 'm', 'LineWidth', 3);
                else
                    set(d.trajInMain, {'XData', 'YData', 'ZData'}, {x,y,z});
                end
            end
        end
        %-----------------------------------------------------------------
        function moveWellXY(d, src, event)
            p = d.trajInMainPnt.Position;
            [x,y] = deal(d.lineOnSlice.XData, d.lineOnSlice.YData);
            d.lineOnSlice.XData = x + (p(1)-x(1));
            d.lineOnSlice.YData = y + (p(2)-y(1));
            d.lineOnSlice.resetLinePoints();
            d.lineOnSlice.updateSlice();
        end
        %-----------------------------------------------------------------
        function layout(d, ~, ~)
            [mw, sp] = deal(d.layoutParams.menuWidth, d.layoutParams.itemSpace);
            fip    = d.Figure.Position;
            mPos   = d.Menu.Position;
            mPos   = [5, fip(4)-mPos(4)-1, mw, mPos(4)];
            aw = floor((fip(3)-mw-4*sp)/2);
            aPos1   = [mw+2*sp, 2*sp, max(10, aw), max(10, fip(4)-3*sp)];
            aPos2   = [mw+3*sp+aw, 2*sp, max(10, aw), max(10, fip(4)-3*sp)];
            d.Menu.Position = mPos;
            d.Axes(1).Position = aPos1;
            d.Axes(2).Position = aPos2;
            if numel(d.Axes) == 3
                mrg   = sp;
                aPos3 = [mrg*[1 1.5], mw-mrg, round(3*(mw-mrg)/4)];
                d.Axes(3).Position = aPos3;
            end
                
        end
        %-----------------------------------------------------------------
        function launchDiagnostics(d, ~, ~)
            n = numel(d.WNew)+1;
            [models, wells, names] = deal(cell(1,n));
            if n>1
                % trim off suffix
                % add closed for first well-list
                newIx = numel(d.W) + 1;%(1:(n-1));
                tmp = [d.W; d.WNew(1)];
                tmp(end).name = 'NEW-WELL';
                [tmp(newIx).status] = deal(false);
                [tmp(newIx).cstatus] = deal(false(size(tmp(end).cstatus)));
                [tmp(newIx).type] = deal('rate');
                [tmp(newIx).val] = deal(0);
                [models{1}, wells{1}]  = deal(d.model, tmp);
                names{1} = 'Base';
                for k = 2:n
                    tmp = [d.W; d.WNew(k-1)];
                    tmp(end).name = 'NEW_WELL';
                    [models{k}, wells{k}]  = deal(d.model, tmp);
                    names{k} = ['case ', num2str(k-1)];
                end
                DiagnosticsViewer(models,wells, 'modelNames', names, 'includeAverage', false);
            end
        end
        %-----------------------------------------------------------------
        function d = switchInteractionOff(d)
            modes = fieldnames(d.interactiveModes);
            for k =1 : numel(modes)
                d.interactiveModes.(modes{k}).enable = 'off';
            end
        end
        %-----------------------------------------------------------------
        function showGradientVectors(d)
            if ~isempty(d.gradientVectors) && any(isvalid(d.gradientVectors))
                delete(d.gradientVectors);
            end
            x = d.lineOnSlice.XData(:);
            y = d.lineOnSlice.YData(:);
            z = d.lineOnSlice.ZData(:);
            g = d.objValues{d.caseNo}.gradient.position{d.wellNo};
            gsc = g*d.objValues{d.caseNo}.gradPlotScale{d.wellNo};
            [u,v,w] = deal(gsc(1:3:end), gsc(2:3:end), gsc(3:3:end));
            cap = @(x,lims)max(lims(1),min(lims(2),x));
            u = cap(x+u, d.Axes(1).XLim)-x;
            v = cap(y+v, d.Axes(1).YLim)-y;
            w = cap(z+w, d.Axes(1).ZLim)-z;
            hold(d.Axes(1), 'on')
            d.gradientVectors(1) = quiver3(d.Axes(1),x,y,z,u,v,w, 'k', 'LineWidth', 2, ...
                'AutoScale', 'off', 'PickableParts', 'none', 'MaxHeadSize', .5);
            ngr = sqrt(sum(reshape(g, 3, []).^2));
            str = arrayfun(@(x)sprintf('  %1.2d $/m', x), ngr, 'UniformOutput', false);
            [d.gradientVectors(2 + (1:numel(x)))] = ...
                text(d.Axes(1), x+u/2,y+v/2,z+w/2, str, 'BackgroundColor', [.92 .92 .92], ...
                    'PickableParts', 'none');
            hold(d.Axes(1), 'off')
            %hold(d.Axes(2), 'on')
            d.gradientVectors(2) = copyobj(d.gradientVectors(1),d.Axes(2));
            %d.gradientVectors(2).Parent = d.Axes(2);
            %hold(d.Axes(2), 'off')
        end
        %-----------------------------------------------------------------
        function moveTrajAlongGradient(d, val)
            if isvalid(d.lineOnSlice)
                gv = d.gradientVectors(1);
                cap = @(x,lims)max(lims(1),min(lims(2),x));
                d.lineOnSlice.XData = cap(gv.XData + val*gv.UData, d.Axes(1).XLim);
                d.lineOnSlice.YData = cap(gv.YData + val*gv.VData, d.Axes(1).YLim);
                d.lineOnSlice.ZData = cap(gv.ZData + val*gv.WData, d.Axes(1).ZLim);
                d.lineOnSlice.resetLinePoints();
                %ngr = sqrt(sum(reshape(d.objValues{d.caseNo}.gradient.position{d.wellNo}, 3, []).^2));
                %dp  = val*[gv.UData(:), gv.VData(:), gv.WData(:)].';
                %dp  = sqrt(sum(dp.^2));
                %inc = sum(ngr.*dp);
                %d.proxySelector.increaseText.String = ...
                %    sprintf('%1.2d $ / %2.1d %%', inc, 100*inc/d.objValues{d.caseNo}.value);
            end
        end
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function traj = getTrajectory(d)
pxy = [d.Line3D.trajectory.XData(:), d.Line3D.trajectory.YData(:)];
p   = [d.Line2D.trajectory.XData(:), d.Line2D.trajectory.YData(:)];
lenxy = [0; sqrt(dot(diff(pxy), diff(pxy), 2))];
%traj = [interp1(fac.cumlength, fac.segments, p(:,1)), p(:,2)];
traj = [interp1(cumsum(lenxy), pxy, p(:,1), 'linear', 'extrap'), p(:,2)];
end
%--------------------------------------------------------------------------
function wp= setSpecialWellPlotProps(wp)
tp = {'injectors', 'producers'};
for k =1:numel(tp)
    % don't show well-names outside region
    set([wp.(tp{k}).label], 'Clipping', 'on');
    % don't interfere with interactivity
    set([wp.(tp{k}).connector], 'PickableParts', 'none');
    set([wp.(tp{k}).body], 'PickableParts', 'none');
    set([wp.(tp{k}).label], 'PickableParts', 'none');
end
% 'shade' color
%wp.injectorColor = 1-.5*(1-wp.injectorColor);
%wp.producerColor = 1-.5*(1-wp.producerColor);
end
%--------------------------------------------------------------------------
function w = updateWellTrajectory(model, w, traj)
if isgraphics(traj, 'Line')
    traj = [traj.XData(:), traj.YData(:), traj.ZData(:)];
end
status = w.status;
opts = {'Name', w.name, 'Type', w.type, 'Val', w.val, 'Sign', w.sign, ...
        'lims', w.lims, 'Radius', w.r(1), 'compi', w.compi};
w  = addWellFromTrajectory([], model.G, model.rock, traj, ...
                    'exteriorFaceCorrection', true, opts{:});
w.status = status;
w.trajectory = traj;
if isempty(w.cells)
    nc = 0;
else
    nc = numel(w.cells);
end
w.cell_origin = ones(nc, 1);
end
%--------------------------------------------------------------------------
function plotObjectiveValues(ax, res, nms)
%figure, ax = gca;
cla(ax, 'reset')
title(ax, 'Proxy values')
v = cellfun(@(r)r.value, res);
plot(ax, v, '--o', 'MarkerSize', 10, 'LineWidth', 2);
ax.XTick = 1:numel(v);
ax.XTickLabel = nms;
set(ax, 'FontSize', 10);
ax.XLim = [.5, numel(v)+.5];
ax.YTickLabelRotation = 45;
ax.XTickLabelRotation = 45;
grid(ax, 'on');
end
%--------------------------------------------------------------------------
function ps = setupPositionControl(G, w, ls)
pSize = [16 16 8];
cPnts = [ls.XData(:), ls.YData(:), ls.ZData(:)];
ps    = WellPositionControl(G, 'w', w, 'perturbationSize', pSize, ...
            'controlPoints', cPnts, 'nPoints', size(cPnts,1));
end
%--------------------------------------------------------------------------
function sc = getGradPlotScaling(g, axLims)
fac = .2;
axLims = vertcat(axLims{:});
axLims = repmat(axLims, numel(g)/3, 1);
rl     = max(abs(g)./diff(axLims, [], 2));
sc     = fac/rl;
end
%--------------------------------------------------------------------------
function data = getCellData(model, state)
tmp.PORO = model.rock.poro;
tmp.PERMX = model.rock.perm(:,1);
tmp.PERMY = model.rock.perm(:,2);
tmp.PERMZ = model.rock.perm(:,3);
tmp.DEPTH = model.G.cells.centroids(:,3);
tmp.PORV  = poreVolume(model.G,model.rock);
if ~isempty(state)
    if model.water
        tmp.SWATINIT = state.s(:, model.getPhaseIndex('W'));
    end
    if model.oil
        tmp.SOILINIT = state.s(:, model.getPhaseIndex('O'));
    end
    if model.gas
        tmp.SGASINIT = state.s(:, model.getPhaseIndex('G'));
    end
end
names = fieldnames(tmp);
vals  = applyFunction(@(nm)tmp.(nm), names);
lims  = applyFunction(@(v)[min(v)*(1-sqrt(eps)), max(v)*(1+sqrt(eps))], vals);
data  = struct('name', names, 'values', vals, 'limits', lims);
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
function plotControlGradients(ax, W, g)
assert(numel(W) == numel(g), 'Non-mathcing gradient')
isp = strcmp({W.type}, 'bhp');
cla(ax, 'reset')
%hold(ax,'on')
%ax.YAxisLocation = 'right';
yyaxis(ax, 'right')
v = g.*isp(:)*barsa;
bar(ax, v);
set(ax, 'YLim', [-1 1]*max(abs(v)));
%xticklabels({W.name});
%xtickangle(45)
ylabel(ax,'Gradients [$/bar]');
yyaxis(ax, 'left')
%ax.YAxisLocation = 'left';
v = g.*(~isp(:))/day;
bar(ax, v);
ax.YLim = [-1 1]*max(abs(v));
xticklabels(ax,{W.name});
ax.XTickLabelRotation = 45;
ylabel(ax,'Gradients [$(m^3/day)^{-1}]')
grid(ax, 'on')
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


 %%-----------------------------------------------------------------
        %         function xyLine(d, ~, ~)
        %             w = d.WNew(d.wellNo);
        %             pnts = w.trajectory(:, 1:2);
        %             pnts = uniquetol(pnts, 'ByRows', true);
        %             if size(pnts, 1) == 1
        %                 dx = d.Axes(1).XLim;
        %                 pnts = [pnts(1)+[dx/5; -dx/5], pnts(2)*ones(2,1)];
        %             end
        %             if ~isempty(d.Line3D)&&isvalid(d.Line3D)
        %                 d.switchInteractionOff();
        %                 delete(d.Line3D);
        %             end
        %             d.Line3D = InteractiveLine('XData', pnts(:,1), 'YData', pnts(:,2), 'Parent', d.Axes(1));
        %             d.Line3D.Callback = @d.crossSection;
        %         end
        
        %%-----------------------------------------------------------------
        %         function putXYLine(d, ~, ~)
        %             xx = d.Axes(1).XLim;
        %             yy = d.Axes(1).YLim;
        %             pnts = [mean(xx) yy(1); mean(xx) yy(2)];
        %             if ~isempty(d.Line3D)&&isvalid(d.Line3D)
        %                 delete(d.Line3D);
        %             end
        %             if ~isempty(d.Line2D)&&isvalid(d.Line2D)
        %                 delete(d.Line2D);
        %             end
        %             d.Line3D = InteractiveLine('XData', pnts(:,1), 'YData', pnts(:,2), 'Parent', d.Axes(1));
        %             d.crossSection();
        %             d.Line3D.Callback = @d.crossSection;
        %             d.empty2DLine();
        %             d.Line2D.smoothing = 'pchip';
        %             d.Axes(2).XLimMode = 'auto';
        %             d.Axes(2).YLimMode = 'auto';
        %         end
        %-----------------------------------------------------------------