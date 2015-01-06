function plotWellSols(wellsols, varargin)

    if mod(numel(varargin), 2) == 1
        timesteps = varargin{1};
        hasTimesteps = true;
        varargin = varargin(2:end);
    else
        hasTimesteps = false;
        timesteps = [];
    end
    
    if isa(wellsols{1}, 'struct')
        % Single input, wrap in cell
        wellsols = {wellsols};
    end
    
    if isa(timesteps, 'double')
        % Single input, wrap in cell
        timesteps = {timesteps};
    end
    
    % Single set of matching timesteps for multiple well sols
    if numel(timesteps) ~= numel(wellsols)
        assert(numel(timesteps) == 1);
        tmp = cell(size(wellsols));
        [tmp{:}] = deal(timesteps{1});
        timesteps = tmp; clear tmp
    end
    
    % Timesteps are either cumulative or individual timesteps. Try to
    % detect if timesteps are actually decreasing or repeated, even though
    % it is not officially supported.
    for ind = 1:numel(timesteps)
        if any(diff(timesteps{ind}) <= 0)
            timesteps{ind} = cumsum(timesteps{ind});
        end
    end
    opt = struct('lowermargin', .1, ...
                 'plotwidth',   .6, ...
                 'linewidth',    2, ...
                 'field',       'bhp', ...
                 'datasetnames', {{}});
    [opt, plotvararg] = merge_options(opt, varargin{:});
    
    % Grab first wellSol and assume it is representative 
    samplews = wellsols{1}{1}(1);
    fn = getNamesFromWS(samplews);
    % Check that default field actually exists
    fnIndex = find(strcmpi(fn, opt.field));
    if isempty(fnIndex)
        warning(['Unknown field ''', opt.field, '''.']);
        fnIndex = 1;
    end
    
    wellnames = arrayfun(@(x) x.name, wellsols{1}{1}, 'UniformOutput', false);
    
    if isempty(opt.datasetnames)
        % If datasets are not actually named, just assign them data1,
        % data2, ..., datan for convenience.
        opt.datasetnames = arrayfun(@(x) ['data', num2str(x)],...
                                   1:numel(wellsols),...
                                   'UniformOutput', false);
    end
    % Make a figure that's wider than the default.
    df = get(0, 'DefaultFigurePosition');
    fh = figure('Position', df.*[1 1 1.75 1]);
    % We want to ensure that the lines are nice and pretty even if the user
    % has messed with the defaults.
    set(fh, 'Renderer', 'painters');
    
    % Options can alter the amount of space the plot itself takes up.
    lm = opt.lowermargin;
    pw = opt.plotwidth;
    % Somewhat magic numbers because the gui in matlab has some magic
    % constants itself.
    plotaxis  = subplot('Position', [.75*lm, lm, pw, 1-2*lm]);
    ctrlpanel = uipanel('Position', [lm+pw, .25*lm, 1-1.25*lm-pw,  1-1.25*lm], ...
                        'Title',  'Selection');
    
    % Left column
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Well property', ...
              'Position',[.01 .9 .45 .1]);
    % Field selection (bhp, water rate etc) 
    fieldsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'popup',...
              'Value', fnIndex, ...
              'String', fn, 'Callback', @drawPlot, ...
              'Position',[.01 .85 .45 .1]);
    % Select between metric and field units
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Unit system', ...
              'Position',[.51 .9 .45 .1]);

    unitsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'popup',...
              'String', {'SI', 'field'}, 'Callback', @drawPlot, ...
              'Position',[.51 .85 .45 .1]);
          
    % Right column
    % Select active wells for plotting
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Well selection', ...
              'Position',[.01 .8 .45 .05]);
    wellsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'listbox', 'Max', 1e9, 'Min', 1,...
              'String', wellnames, 'Callback', @drawPlot, ...
              'Position',[.01 .2 .45 .6]);
    
    % Select various options
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Misc', ...
              'Position',[.51 .8 .45 .05]);
          
    bg = uibuttongroup('Units', 'normalized', 'Parent', ctrlpanel,...
              'Position',[.51 .2 .45 .6]);

    % Show minor grid to make plot easier to read.
    gridtog = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', true, ...
              'String','Grid on', 'Callback', @drawPlot, ...
              'Position',[.01 .9 .95 .1]);
    % Log transform x axis
    logx = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', ...
              'String','Log (x)', 'Callback', @drawPlot, ...
              'Position',[.01 .8 .95 .1]);
    % Log transform y axis
    logy = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', ...
              'String','Log (y)', 'Callback', @drawPlot, ...
              'Position',[.01 .7 .95 .1]);
    % Draw markers
    hasmarker = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', ...
              'String','Markers', 'Callback', @drawPlot, ...
              'Position',[.01 .6 .95 .1]);
    % Show legend with well and dataset names
    useleg = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', 1, ...
              'String','Legend', 'Callback', @drawPlot, ...
              'Position',[.01 .5 .95 .1]);
    legh = nan;
    % Cumulative sum - used for for example water production, but does not
    % necessarily make sense for pressure.
    csum = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', 0, ...
              'String','Cumulative sum', 'Callback', @drawPlot, ...
              'Position',[.01 .4 .95 .1]);
    % Take abs value of data.
    abst = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', 0, ...
              'String','Absolute value', 'Callback', @drawPlot, ...
              'Position',[.01 .3 .95 .1]);
          
    if hasTimesteps
        % Toggle to use timesteps for spacing, otherwise the x nodes will
        % be equidistant.
        showdt = uicontrol('Units', 'normalized', 'Parent', bg,...
                  'Style', 'checkbox', 'Value', hasTimesteps,...
                  'String','Use timesteps', 'Callback', @drawPlot, ...
                  'Position',[.01 .2 .95 .1]);
    end
    % Line width of the plot
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Plot line width', ...
              'Position',[.01 .14 .99 .05]);
          
    wsl = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
                'Style', 'slider', 'Min', 0, 'Max', 10,...
                'Value', opt.linewidth, 'Callback', @drawPlot, ...
                'Position',[.01 .1 .99 .05]);
    % Figure -> jpg/png/eps
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
                'String','Save as image', 'Callback', @savePlotAsFigure, ...
                'Position',[.01 .01 .485 .08]);
    % Save the plotted data as variables in the workspace.
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
                'String','Save to workspace', ...
                'Callback', @saveDataToWorkSpace, ...
                'Position',[.51 .01 .485 .08]);

    % Initialize multiscope variables
    prevWells = {};
    wells = {};
    fld = '';
    currentdata = {};
    % Initial plot sets up these fields for us and puts something on the screen.
    drawPlot([], []);
    
    function drawPlot(src, event, varargin)
        fld = getFieldString(fieldsel, true);
        wells = getFieldString(wellsel, false);
        
        ndata = numel(wellsols);
        nw = numel(wells);
        
        if nw == 1
            ncolors = ndata;
        else
            ncolors = nw;
        end
        
        if ncolors < 8
            % Relatively few wells, use matlab default
            cmap = lines(ncolors);
        else
            % Use colorcube for more lines
            cmap = colorcube(ncolors+2);
        end
        linestyles = {'--', ':', '-', '-.'};
        
        if get(hasmarker, 'Value')
            markerstyles = {'o', '.', 'd', '*'};
        else
            markerstyles = {''};
        end
            
        
        axis(plotaxis);
        cla; hold on;
        
        l = {};

        tit = '';
        currentdata = cell(ndata, 1);
        for i = 1:ndata
            for j = 1:nw
                wname = wells{j};
                line = linestyles{mod(i-1, numel(linestyles)) + 1};
                m = markerstyles{mod(i-1, numel(markerstyles)) + 1};
                
                d = getData(wname, wellnames, fld, wellsols{i});
                if hasTimesteps && get(showdt, 'Value')
                    x = timesteps{i}/day;
                    xlabel('Time (days)')
                    xunit = day;
                else
                    x = 1:numel(d);
                    xunit = 1;
                    xlabel('Step #')
                end
                
                if get(csum, 'Value')
                    d = cumtrapz(x*xunit, d);
                end
                
                if get(abst, 'Value')
                    d = abs(d);
                end
                
                [tit, d, yl] = getWellUnit(d, fld, getFieldString(unitsel, true));
                ylabel(yl);

                linew = get(wsl, 'Value');
                if linew == 0;
                    % Draw no lines
                    line = '';
                    linew = 1;
                end
                
                if isfield(wellsols{i}{1}, 'status')
                    % Mask away inactive data points
                    status = getData(wname, wellnames, 'status', wellsols{i});
                    d(status == 0, :) = nan;
                end
                if nw == 1 
                    c = cmap(i, :);
                else
                    c = cmap(j, :);
                end
                plot(x, d, [m, line], 'LineWidth', linew, 'color', c, plotvararg{:});
                
                tmp = wname;
                if ndata > 1
                    tmp = [tmp, ' (' opt.datasetnames{i}, ')']; %#ok
                end
                l = [l; tmp];
                currentdata{i} = [currentdata{i}, d];
            end
        end
        title(tit)
        
        if ~isempty(legh) && ishandle(legh) && numel(wells) == numel(prevWells)
            lpos = get(legh, 'Position');
        else
            lpos = 'NorthEast';
        end
        
        if get(useleg, 'Value')
            legh = legend(l, 'Location', lpos);
        elseif ishandle(legh)
            delete(legh);
            legh = nan;
        end
        if get(gridtog, 'Value')
            grid on
        end
        
        if get(logy, 'Value')
            set(plotaxis, 'YScale', 'log')
        else
            set(plotaxis, 'YScale', 'linear')
        end
        
        if get(logx, 'Value')
            set(plotaxis, 'XScale', 'log')
        else
            set(plotaxis, 'XScale', 'linear')
        end
        
        axis tight
        prevWells = wells;
    end

    function savePlotAsFigure(src, event)
        plotname = fld;
        if numel(wells) > 0
            plotname = [plotname, sprintf('_%s', wells{:})];
        end
        ext = '*.jpg; *.tif; *.png; *.gif; *.eps; *.pdf';
        [file, path] = ...
        uiputfile({ext, ['Images (', ext ,')']; ...
                  '*.*', 'All Files'} ,...
                  'Save Image', [plotname, '.jpg']);
       
       if isnumeric(file) && file == 0
           % User aborted action
           return
       end
       % Need this to make a copy of the axis
       dap = get(0, 'DefaultAxesPosition');
       
       dims = get(fh, 'Position');
       newdims = dims.*[1, 1, pw + dap(1), 1];
       tmpfig = figure('Position', newdims);
       ax = copyobj(plotaxis, tmpfig);
       % Approx same dimensions as the gui reference
       set(ax, 'Position', dap)
       
       if exist('export_fig', 'file') == 2
           % Export fig is installed and we can call it
           export_fig([path, file], '-transparent');
           close(tmpfig);
           return
       end
       
       [fileext, fileext, fileext] = fileparts(file); %#ok
       fileext = fileext(2:end);
       
       switch(lower(fileext))
           case {'jpg', 'jpeg'}
               popt = '-djpeg';
           case {'tif', 'tiff'}
               popt = '-dtiff';
           case 'png'
               popt = '-dpng';
           case 'pdf'
               popt = '-dpdf';
           case 'eps'
               popt = '-depsc';
           otherwise
               error(['Unable to save file. Unknown file extension ', fileext]); 
               close(tmpfig);
       end
       print(tmpfig, popt, [path, file]);
       close(tmpfig);
    end

    function saveDataToWorkSpace(src, event)
        clean = @(x) x(~isspace(x));
        varnames = cellfun(@(x) clean(x), opt.datasetnames, 'UniformOutput', false);
        
        % All current datasets + names of selected wells
        export2wsdlg([opt.datasetnames, 'Selected well names'], ...
                     [varnames, 'wellnames'], ...
                     {currentdata{:}, wells})
    end
end

function [tit, d, yl] = getWellUnit(d, fld, usys)
    isMetric = strcmpi(usys, 'si');
    
    yl = '';
    tit = fld;
    switch lower(fld)
        case {'qws', 'qos', 'qgs', 'rate', 'qts', 'qwr', 'qgr', 'qor', 'qtr'}
            switch lower(fld)
                case {'qos', 'qor'}
                    ph = 'oil';
                case {'qws', 'qwr'}
                    ph = 'water';
                case {'qgs', 'qgr'}
                    ph = 'gas';
                otherwise
                    ph = 'total';
            end
            
            if numel(fld) == 3 && lower(fld(3)) == 'r'
                tmp = 'reservoir';
            else
                tmp = 'surface';
            end
            tit = [fld, ': Well ', tmp, ' rate (', ph, ')'];
            if isMetric
                yl = 'm^3/s';
            else
                yl = 'stb/day';
                d = convertTo(d, stb/day);
            end
        case {'bhp', 'pressure'}
            tit = [fld, ': Bottom hole pressure'];
            if ~isMetric
                yl = 'Pascal';
            else
                yl = 'Barsa';
                d = convertTo(d, barsa);
            end
        case 'gor'
            tit = [fld, ': Gas/oil ratio at surface conditions'];
        case {'ocut', 'wcut', 'gcut'}
            switch lower(fld(1))
                case 'o'
                    t = 'Oil';
                case 'g'
                    t = 'Gas';
                case 'w'
                    t = 'Water';
            end
            tit = [fld, ': ', t, ' fraction at reservoir conditions'];
        case 'sign'
            tit = [fld, ': Well sign (+1 for producer, -1 for injector'];
        case 'val'
            tit = [fld, ': Well control value'];
        case 'cdp'
            tit = [fld, ': Pressure drop from reference depth to first perforation'];
            if ~isMetric
                yl = 'Pascal';
            else
                yl = 'Barsa';
                d = convertTo(d, barsa);
            end
        otherwise
            disp('Unknown well field - no unit found');
    end
end

function fn = getNamesFromWS(ws)
    fn = fieldnames(ws);
    % Filter non-numerical values and values with multiple values per well
    fn = fn(cellfun(@(x) isnumeric(ws.(x)) && numel(ws.(x)) == 1, fn));
    fn = fn(end:-1:1);
end

function f = getFieldString(handle, issingle)
    if nargin == 1
        issingle = false;
    end
    s = get(handle, 'String');
    f = s(get(handle, 'Value'));
    if issingle
        f = f{1};
    end
end

function d = getData(wellname, wellnames, field, ws)
    ind = find(strcmpi(wellname, wellnames));
    assert(numel(ind) == 1, 'Multiple wells with same name!');
    d = cellfun(@(x) x(ind).(field), ws);
end
