function varargout = plotWellSols(wellsols, varargin)
%Plot well solutions from AD-solvers
%
% SYNOPSIS:
%   plotWellSols(wellSols, time);
%   plotWellSols(wellSols);
%
% DESCRIPTION:
%   Open interactive plotting interface for well solutions.
%
% REQUIRED PARAMETERS:
%   wellSols - Cell array of `nstep` by 1, each containing a uniform struct
%              array of well solution structures. For example, the first
%              output from `simulateScheduleAD`. Can also be a cell array 
%              of such cell arrays, for comparing multiple simulation
%              scenarios.
%
%   time     - (OPTIONAL) The time for each timestep. If not provided, the
%              plotter will use step number as the x axis intead. If
%              wellSols is a cell array of multiple datasets, time should
%              also be a cell array, provided not all datasets use the same
%              timesteps.
%
% OPTIONAL PARAMETERS:
%   'field'        -  Initial field for plotting (default: 'bhp').
%
%   'linestyles'   - Cell array of line styles used for different datasets.
%
%   'markerstyles' - Marker array of line styles used for different
%                    datasets. 
%
%   'datasetnames' - A cell array of dataset names used for the legend when
%                    plotting multiple datasets.
%
%   'timescale'    - A string for the default choice for axis time-scale. A
%                    string which matches either choice:
%                    'days', 'minutes', 'seconds', 'hours', 'years'
%
%   toggleOn       - A string followed by true/false to enable/disable the
%                    toggles in the Misc panel. The strings are: 'grid',
%                    'logx', 'logy', 'marker', 'legend', 'cumsum', 'abs',
%                    'zoom', 'stairs', and 'ctrl'
%
% RETURNS:
%   fh     - figure handle to plotting panel
%
%   inject - function handle used to dynamically inject new datasets into
%            the viewer (for example, from a running simulation). Same
%            syntax as the base function, but does not support additional
%            varargin.
%
% SEE ALSO:
%   `simulateScheduleAD`

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
    opt = struct('lowermargin',   .1, ...
                 'plotwidth',     .6, ...
                 'linewidth',     2, ...
                 'field',         'bhp', ...
                 'SelectedWells', 1, ...
                 'linestyles',    {{'-', '--', '-.', ':'}}, ...
                 'markerstyles',  {{'o', 'x', 'd', 's', 'd', '^', 'v', '>', '<', '*'}}, ...
                 'timescale',     'days', ...
                 'grid',          true,   ...
                 'logx',          false,  ...
                 'logy',          false,  ...
                 'marker',        false,  ...
                 'legend',        true,   ...
                 'cumsum',        false,  ...
                 'abs',           true,   ...
                 'zoom',          false,  ...
                 'stairs',        false,  ...
                 'ctrl',          false,  ...
                 'figure',        [],     ...
                 'datasetnames', {{}});
    if (iscell(wellsols) && isfield(wellsols{1}, 'SimulatorSetup')) ...
            || isfield(wellsols, 'SimulatorSetup')
        % We were given one or more packed problems and we get everything
        % from those
        problems = wellsols;
        if isstruct(problems)
            problems = {problems};
        end
        [wellsols, ~, ~, names, timesteps] = getMultiplePackedSimulatorOutputs(problems,...
                    'readFromDisk', false, 'readWellSolsFromDisk', true);
        hasTimesteps = true;
        opt.datasetnames = names;
    else
        % We were given individual arguments, figure out what is what.
        if mod(numel(varargin), 2) == 1
            timesteps = varargin{1};
            hasTimesteps = true;
            varargin = varargin(2:end);
        else
            hasTimesteps = false;
            timesteps = [];
        end
    end

    if isa(wellsols{1}, 'struct')
        % Single input, wrap in cell
        wellsols = {wellsols};
    end
    ndata = numel(wellsols);
    warning('ON', 'mrst:badcumsum');
    timesteps = validateTimesteps(wellsols, timesteps);
    
    % Timesteps are either cumulative or individual timesteps. Try to
    % detect if timesteps are actually decreasing or repeated, even though
    % it is not officially supported.
    for ind = 1:numel(timesteps)
        if any(diff(timesteps{ind}) <= 0)
            timesteps{ind} = cumsum(timesteps{ind});
        end
    end
    [opt, plotvararg] = merge_options(opt, varargin{:});
    
    % Grab first wellSol and assume it is representative
    if isempty(wellsols{1}{1})
        disp('No wellSols present in input.');
        return
    else
        samplews = wellsols{1}{1}(1);
    end
    fn = getNamesFromWS(samplews);
    % Check that default field actually exists
    fnIndex = find(strcmpi(fn, opt.field));
    if isempty(fnIndex)
        warning(['Unknown field ''', opt.field, '''.']);
        fnIndex = 1;
    end
    all_well_names = cell(numel(wellsols), 1);
    for ind = 1:numel(wellsols)
        all_well_names{ind} = arrayfun(@(x) x.name, wellsols{ind}{1}, 'UniformOutput', false);
    end
    
    wellnames = all_well_names{1};
    
    if isempty(opt.datasetnames)
        % If datasets are not actually named, just assign them data1,
        % data2, ..., datan for convenience.
        opt.datasetnames = arrayfun(@(x) ['data', num2str(x)],...
                                   1:numel(wellsols),...
                                   'UniformOutput', false);
    end
    % Make a figure that's wider than the default.
    df = get(0, 'DefaultFigurePosition');
    if isempty(opt.figure) || ~ishandle(opt.figure)
        fh = figure('Position', df.*[1 1 1.75 1]);
    else
        fh = opt.figure;
    end
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
    xmargin = 0.01;
    midmargin = xmargin/2;
    blocksz = (1 - 2*xmargin - midmargin)/2;
    leftOffset = blocksz + xmargin + midmargin;
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Well property', ...
              'Position',[xmargin, .9, blocksz, .1]);
    % Field selection (bhp, water rate etc) 
    fieldsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'popup',...
              'Value', fnIndex, ...
              'String', fn, 'Callback', @drawPlot, ...
              'Position',[xmargin, .85, blocksz, .1]);
    % Select between metric and field units
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Unit system', ...
              'Position',[leftOffset, .9, blocksz, .1]);

    unitsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'popup',...
              'String', {'Metric', 'SI', 'Field', 'Lab'}, 'Callback', @drawPlot, ...
              'Position',[leftOffset, .85, blocksz/2, .1]);
    timechoices = {'Years', 'Days', 'Hours', 'Minutes', 'Seconds'};
    timescales = [year(), day(), hour(), minute(), second()];
    
    timesel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'popup',...
              'Value', find(strcmpi(opt.timescale, timechoices)), ...
              'String', timechoices, ...
              'Callback', @drawPlot, ...
              'Position',[leftOffset + blocksz/2, .85, blocksz/2, .1]);
    % Right column
    % Select active wells for plotting
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Well selection', ...
              'Position',[xmargin, .8, blocksz, .05]);
  
      if ndata > 1
          wellsel_size = [xmargin, .5, blocksz, .3];
          uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Dataset selection', ...
              'Position',[xmargin, .45, blocksz, .05]);
          datasel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'listbox', 'Max', 1e9, 'Min', 1,...
              'String', opt.datasetnames, 'Callback', @drawPlot, ...
              'Value', 1:numel(opt.datasetnames), ...
              'Position', [xmargin, .2, blocksz, .25]);
      else
          wellsel_size = [xmargin, .2, blocksz, .6];
      end
    wellsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'listbox', 'Max', 1e9, 'Min', 1,...
              'Value', opt.SelectedWells, ...
              'String', wellnames, 'Callback', @drawPlot, ...
              'Position', wellsel_size);
    
    % Select various options
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Misc', ...
              'Position',[leftOffset, .8, blocksz, .05]);
          
    bg = uibuttongroup('Units', 'normalized', 'Parent', ctrlpanel,...
              'Position',[leftOffset, .2, blocksz, .6]);

    % Show minor grid to make plot easier to read.
    toggle_h = 1./(10 + double(hasTimesteps || nargout > 1));
    gridtog = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', opt.grid, ...
              'String','Grid on', 'Callback', @drawPlot, ...
              'Position',[xmargin, 1 - toggle_h, 1-xmargin .1]);
    % Log transform x axis
    logx = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'value', opt.logx, ...
              'String','Log (x)', 'Callback', @drawPlot, ...
              'Position',[xmargin, 1 - 2*toggle_h, 1-xmargin .1]);
    % Log transform y axis
    logy = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'value', opt.logy, ...
              'String','Log (y)', 'Callback', @drawPlot, ...
              'Position',[xmargin, 1 - 3*toggle_h, 1-xmargin .1]);
    % Draw markers
    hasmarker = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox',  'value', opt.marker, ...
              'String','Markers', 'Callback', @drawPlot, ...
              'Position',[xmargin, 1 - 4*toggle_h, 1-xmargin .1]);
    % Show legend with well and dataset names
    useleg = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', opt.legend, ...
              'String','Legend', 'Callback', @drawPlot, ...
              'Position',[xmargin, 1 - 5*toggle_h, 1-xmargin .1]);
    legh = nan;
    % Cumulative sum - used for for example water production, but does not
    % necessarily make sense for pressure.
    csum = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', opt.cumsum, ...
              'String','Cumulative sum', 'Callback', @drawPlot, ...
              'Position',[xmargin, 1 - 6*toggle_h, 1-xmargin .1]);
    % Take abs value of data.
    abst = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', opt.abs, ...
              'String','Absolute value', 'Callback', @drawPlot, ...
              'Position',[xmargin, 1 - 7*toggle_h, 1-xmargin .1]);
          
    % Zoom to data range
    zoomt = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', opt.zoom, ...
              'String','Zoom to data', 'Callback', @drawPlot, ...
              'Position',[xmargin, 1 - 8*toggle_h, 1-xmargin .1]);
    stairplot = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', opt.stairs,...
              'String','Stair-step plot', 'Callback', @drawPlot, ...
              'Position',[xmargin, 1 - 9*toggle_h, 1-xmargin .1]);
    ctrlplot = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', opt.ctrl,...
              'String','Show control changes', 'Callback', @drawPlot, ...
              'Position',[xmargin, 1 - 10*toggle_h, 1-xmargin .1]);
    if hasTimesteps || nargout > 1
        % Toggle to use timesteps for spacing, otherwise the x nodes will
        % be equidistant.
        showdt = uicontrol('Units', 'normalized', 'Parent', bg,...
                  'Style', 'checkbox', 'Value', hasTimesteps ,...
                  'String','Use timesteps', 'Callback', @drawPlot, ...
                  'Position',[xmargin, 1 - 11*toggle_h, 1-xmargin .1]);
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
                'ButtonDownFcn', @savePlotAsFigure, ...
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
    
    function drawPlot(src, event, varargin) %#ok<*INUSD>
        set(0, 'CurrentFigure', fh);
        fld = getFieldString(fieldsel, true);
        wells = getFieldString(wellsel, false);

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
        linestyles = opt.linestyles;
        
        if get(hasmarker, 'Value')
            markerstyles = opt.markerstyles;
        else
            markerstyles = {''};
        end
            
        
        axis(plotaxis);
        cla; hold on;
        
        l = {};

        tit = '';
        currentdata = cell(ndata, 1);
        Mv = -inf;
        mv =  inf;
        if ndata == 1
            data_iter = 1;
        else
            data_iter = get(datasel, 'Value');
        end
        
        for i = data_iter
            for j = 1:nw
                wname = wells{j};
                line = linestyles{mod(i-1, numel(linestyles)) + 1};
                m = markerstyles{mod(i-1, numel(markerstyles)) + 1};
                
                [d, wpos] = getData(wname, all_well_names{i}, fld, wellsols{i});
                if hasTimesteps && get(showdt, 'Value')
                    timescaleix = get(timesel, 'Value');
                    nowTime = timechoices{timescaleix};
                    timescale = timescales(timescaleix);
                    x = timesteps{i}/timescale;
                    xlabel(['Time [', nowTime, ']'])
                    xunit = day;
                else
                    x = 1:numel(d);
                    xunit = 1;
                    xlabel('Step #')
                end
                
                [tit, d, yl, doCsum] ...
                    = getWellUnit(d, fld, ...
                                  getFieldString(unitsel, true), ...
                                  get(csum,'Value'), ...
                                  hasTimesteps);
                ylabel(yl);
                
                if numel(wellsols{i}) > 0 && isfield(wellsols{i}{1}, 'status')
                    % Mask away inactive data points
                    status = getData(wname, all_well_names{i}, 'status', wellsols{i});
                    active_step = status > 0;
                else
                    active_step = true(size(d));
                end
                
                if doCsum
                    d = cumtrapz(x*xunit, active_step.*d);
                end
                d(~active_step, :) = nan;

                if get(abst, 'Value')
                    d = abs(d);
                end

                linew = get(wsl, 'Value');
                if linew == 0
                    % Draw no lines
                    line = '';
                    linew = 1;
                end

                if nw == 1 
                    c = cmap(i, :);
                else
                    c = cmap(j, :);
                end
                if get(stairplot, 'Value')
                    pltfn = @stairs;
                else
                    pltfn = @plot;
                end
                
                pltfn(x, d, [m, line], 'LineWidth', linew, 'color', c, plotvararg{:});
                Mv = max(Mv, max(d(:)));
                mv = min(mv, min(d(:)));
                if get(ctrlplot, 'Value') && isfield(wellsols{i}{1}, 'type')
                    type = cellfun(@(ws) ws(wpos).type, wellsols{i}, 'UniformOutput', false);
                    [types, tmp, ctrl_id] = unique(type);                  %#ok
                    switched = [1; find(diff(ctrl_id, 1) ~= 0)];
                    for switch_step = 1:numel(switched)
                        six = switched(switch_step);
                        if switch_step == 1
                            sl = type{six};
                        else
                            sl = [type{six}, '->', type{six+1}];
                        end
                        text(x(six), d(six), sl, ...
                            'HorizontalAlignment', 'center', 'BackgroundColor', [1, 1, 1]*0.7,...
                            'LineStyle', line, 'edgecolor', c, 'linewidth', 2);
                    end
                end
                
                tmp = wname;
                if ndata > 1
                    if nw == 1
                        tmp = opt.datasetnames{i};
                    else
                        tmp = [tmp, ' (' opt.datasetnames{i}, ')']; %#ok
                    end
                end
                l = [l; tmp];
                currentdata{i} = [currentdata{i}, d];
            end
        end
        title(tit)
        
        if ~isempty(legh) && ishandle(legh) && numel(wells) == numel(prevWells)
            lpos = get(legh, 'Position');
        else
            % lpos = 'NorthEast';
            lpos = 'Best';
        end
        
        if get(useleg, 'Value')
            legh = legend(l, 'Location', lpos, 'Interpreter', 'none');
        elseif ishandle(legh)
            delete(legh);
            legh = nan;
        end
        if get(gridtog, 'Value')
            grid on
        else
            grid off
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
        if ~get(zoomt, 'Value')
            % We should ensure that zero axis is included
            sM = sign(Mv); sm = sign(mv);
            if sM == sm
                if sM > 0
                    ylm = [0, Mv];
                else
                    ylm = [mv, 0];
                end
            else
                ylm = [min(mv, 0), max(Mv, 0)];
            end
            ylim(ylm + [-eps, eps]);
        end

        prevWells = wells;
    end

    function savePlotAsFigure(src, event) %#ok<INUSL>
        if isempty(event)
            doSave = true;
        else
            doSave = strcmp(event.EventName, 'Action');
        end
        
        if doSave
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
        end
       % Need this to make a copy of the axis
       dap = get(0, 'DefaultAxesPosition');
       dims = get(fh, 'Position');
       newdims = dims.*[1, 1, pw + 2*.75*lm, 1];
       tmpfig = figure('Position', newdims);
       set(tmpfig,'PaperPositionMode','auto');
       if isnan(double(legh))
           ax = copyobj(plotaxis, tmpfig);
       else
           ax = copyobj([legh, plotaxis], tmpfig);
           ax = ax(2);
       end
       % Refresh new axes
       drawnow
       % Approx same dimensions as the gui reference
       set(ax, 'Position', dap)
       if doSave
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
                   close(tmpfig);
                   error(['Unable to save file. Unknown file extension ', fileext]); 
           end
           print(tmpfig, popt, [path, file]);
           close(tmpfig);
       end
    end

    function saveDataToWorkSpace(src, event)
        clean = @(x) x(~isspace(x));
        varnames = cellfun(@(x) clean(x), opt.datasetnames, 'UniformOutput', false);
        
        % All current datasets + names of selected wells
        export2wsdlg([opt.datasetnames, 'Selected well names'], ...
                     [varnames, 'wellnames'], ...
                     {currentdata{:}, wells})
    end

    function injectDataset(ws, steps)
        if nargin == 1
            hasTimesteps = false;
            steps = [];
        else
            hasTimesteps = true;
        end
        wellsols = ws;
        timesteps = validateTimesteps(wellsols, steps);
        drawPlot([], []);
    end
    inject = @(ws, varargin) injectDataset(ws, varargin{:});
    if nargout > 0
        varargout{1} = fh;
        if nargout > 1
            varargout{2} = inject;
        end
    end
end

function [tit, d, yl, doCsum] = getWellUnit(d, fld, usys, isCsum, hasTimesteps)    
    doCsum = isCsum;
    yl = '';
    tit = fld;
    
    if hasTimesteps
        switch lower(usys)
            case 'metric'
                t_unt = day();
                t_str = 'day';
            case 'field'
                t_unt = day();
                t_str = 'day';
            case 'si'
                t_unt = second();
                t_str = 's';
            case 'lab'
                t_unt = hour();
                t_str = 'hour';
        end
    else
        t_unt = 1;
        t_str = 's';
        if isCsum
            warning('mrst:badcumsum', ['Timesteps not provided to plotWellSols,', ...
                    ' cumulative sum will not be accurate']);
            warning('OFF', 'mrst:badcumsum');
        end
    end
    
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
                isRes = true;
                tmp = 'reservoir';
            else
                isRes = false;
                tmp = 'surface';
            end
            if isCsum
                tit = [fld, ': Cumulative ', tmp, ' production (', ph, ')'];
            else
                tit = [fld, ': Well ', tmp, ' rate (', ph, ')'];
            end
            switch lower(usys)
                case 'metric'
                    y_str = 'm^3';
                    y_unit = meter^3;
                case 'field'
                    if isRes
                        y_str = 'stb';
                        y_unit = stb;
                    else
                        if strcmpi(ph, 'gas')
                            y_str = 'scf';
                            y_unit = ft^3;
                        else
                            y_str = 'stb';
                            y_unit = stb;
                        end
                    end
                case 'si'
                    y_str = 'm^3';
                    y_unit = meter^3;
                case 'lab'
                    y_str = 'cm^3';
                    y_unit = (centi*meter)^3;
            end
            if isCsum
                d = convertTo(d, y_unit);
                yl = y_str;
            else
                d = convertTo(d, y_unit/t_unt);
                yl = [y_str, '/', t_str];
            end
            doCsum = isCsum;
        case {'bhp', 'pressure'}
            tit = [fld, ': Bottom hole pressure'];
            [y_unit, yl] = getPressureUnit(usys);
            d = convertTo(d, y_unit);
        case 'gor'
            tit = [fld, ': Gas/oil ratio at surface conditions'];
            doCsum = isCsum;

            switch lower(usys)
               case {'metric', 'si'}
                  yl = 'Sm^3/Sm^3';

               case 'field'
                  yl = 'MScf/stb';
                  d  = convertTo(d, 1000*ft^3 / stb);

               case 'lab'
                  yl = 'Scm^3/Scm^3';
            end
        case 't'
            tit = [fld, ': Temperature'];
            yl  = 'kelvin';
        case {'ocut', 'wcut', 'gcut'}
            switch lower(fld(1))
                case 'o'
                    t = 'Oil';
                case 'g'
                    t = 'Gas';
                case 'w'
                    t = 'Water';
            end
            tit = [fld, ': ', t, ' volume fraction at reservoir conditions'];
        case 'sign'
            tit = [fld, ': Well sign (+1 for injector, -1 for producer)'];
        case 'val'
            tit = [fld, ': Well control value'];
        case 'cdp'
            tit = [fld, ': Pressure drop from reference depth to first perforation'];
            [y_unit, yl] = getPressureUnit(usys);
            d = convertTo(d, y_unit);
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

function [d, ind] = getData(wellname, wellnames, field, ws)
    ind = find(strcmpi(wellname, wellnames));
    assert(numel(ind) == 1, 'Multiple wells with same name!');
    d = cellfun(@(x) x(ind).(field), ws);
end

function timesteps = validateTimesteps(wellsols, timesteps)
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
end

function [unit, str] = getPressureUnit(usys)
    switch lower(usys)
        case 'metric'
            str = 'barsa';
            unit = barsa;
        case 'field'
            str = 'psia';
            unit = psia;
        case 'si'
            str = 'Pascal';
            unit = Pascal;
        case 'lab'
            str = 'atm';
            unit = atm;
    end
end
