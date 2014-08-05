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
        timesteps = {timesteps};
    end
    
    % Grab first and best element for testing
    samplews = wellsols{1}{1}(1);
    fn = getNamesFromWS(samplews);

    wellnames = arrayfun(@(x) x.name, wellsols{1}{1}, 'UniformOutput', false);
    
    opt = struct('lowermargin', .1, ...
                 'plotwidth',   .6, ...
                 'linewidth',    2);
    opt = merge_options(opt, varargin{:});
    
    df = get(0, 'DefaultFigurePosition');
    
    fh = figure('Position', df.*[1 1 1.75 1]);
    
    set(fh, 'Renderer', 'painters');
    
    lm = opt.lowermargin;
    pw = opt.plotwidth;
    
    plotaxis  = subplot('Position', [.5*lm, lm, pw, 1-2*lm]);
    ctrlpanel = uipanel('Position', [lm+pw, lm, 1-1.25*lm-pw,  1-2*lm], ...
                        'Title',  'Selection');
    
    % Top row

    
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Well property', ...
              'Position',[.01 .9 .45 .1]);
          
    fieldsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'popup',...
              'String',fn, 'Callback', @drawPlot, ...
              'Position',[.01 .85 .45 .1]);
          
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Unit system', ...
              'Position',[.51 .9 .45 .1]);
          
    unitsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'popup',...
              'String', {'SI', 'field'}, 'Callback', @drawPlot, ...
              'Position',[.51 .85 .45 .1]);
          
    % Second row
    % Select wells
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Well selection', ...
              'Position',[.01 .8 .45 .05]);
          
    wellsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'listbox', 'Max', 1e9,...
              'String', wellnames, 'Callback', @drawPlot, ...
              'Position',[.01 .2 .45 .6]);
    
    % Select various options
    uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Misc', ...
              'Position',[.51 .8 .45 .05]);
          
    bg = uibuttongroup('Units', 'normalized', 'Parent', ctrlpanel,...
              'Position',[.51 .2 .45 .6]);
    
    gridtog = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', true, ...
              'String','Grid on', 'Callback', @drawPlot, ...
              'Position',[.01 .9 .95 .1]);
          
    logx = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', ...
              'String','Log (x)', 'Callback', @drawPlot, ...
              'Position',[.01 .8 .95 .1]);
          
    logy = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', ...
              'String','Log (y)', 'Callback', @drawPlot, ...
              'Position',[.01 .7 .95 .1]);
          
    hasmarker = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', ...
              'String','Markers', 'Callback', @drawPlot, ...
              'Position',[.01 .6 .95 .1]);
          
    useleg = uicontrol('Units', 'normalized', 'Parent', bg,...
              'Style', 'checkbox', 'Value', 1, ...
              'String','Legend', 'Callback', @drawPlot, ...
              'Position',[.01 .5 .95 .1]);
    legh = nan;
    
    if hasTimesteps
        showdt = uicontrol('Units', 'normalized', 'Parent', bg,...
                  'Style', 'checkbox', 'Value', hasTimesteps,...
                  'String','Use timesteps', 'Callback', @drawPlot, ...
                  'Position',[.01 .4 .95 .1]);
    end
        uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'text',...
              'String', 'Linewidth', ...
              'Position',[.01 .12 .99 .05]);
          
    wsl = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
                'Style', 'slider', 'Min', sqrt(eps), 'Max', 10,...
                'Value', opt.linewidth, ...
                'String','Apply', 'Callback', @drawPlot, ...
                'Position',[.01 .01 .99 .1]);
    
%     % Buttons
%     uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
%               'String','Apply', 'Callback', @drawPlot, ...
%               'Position',[.01 .01 .4 .1]);
%           
%     uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
%               'Style', 'edit', ...
%               'String','--o', 'Callback', @drawPlot, ...
%               'Position',[.51 .01 .4 .1]);
          
    function drawPlot(src, event, varargin)
        fld = getFieldString(fieldsel, true);
        wells = getFieldString(wellsel, false);
        
        ndata = numel(wellsols);
        nw = numel(wells);
        
        cmap = lines(nw);
        linestyles = {'-', '--', ':', '-.'};
        if get(hasmarker, 'Value')
            m = 'o';
        else
            m = '';
        end
            
        
        axis(plotaxis);
        cla; hold on;
        
        l = {};
        for i = 1:ndata
            for j = 1:nw
                wname = wells{j};
                line = linestyles{mod(i-1, numel(linestyles) + 1) + 1};
                
                
                d = getData(wname, wellnames, fld, wellsols{i});
                [d, yl] = getWellUnit(d, fld, getFieldString(unitsel, true));
                ylabel(yl);
                if hasTimesteps && get(showdt, 'Value')
                    x = timesteps{i}/day;
                    xlabel('Time (days)')
                else
                    x = 1:numel(d);
                    xlabel('Step #')
                end
                
                plot(x, d, [m, line], 'LineWidth', get(wsl, 'Value'), 'color', cmap(j, :));
                
                tmp = wname;
                if ndata > 1
                    tmp = [tmp, ' dataset ', num2str(i)]; %#ok
                end
                l = [l; tmp];
            end
        end
        title(fld)
        
        if get(useleg, 'Value')
            legh = legend(l);
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
    end
end

function [d, yl] = getWellUnit(d, fld, usys)
    isMetric = strcmpi(usys, 'si');
    
    yl = '';
    switch lower(fld)
        case {'qws', 'qos', 'qgs', 'rate'}
            if isMetric
                yl = 'Well surface rate (m^3/s)';
            else
                yl = 'Well surface rate (stb/day)';
                d = convertTo(d, stb/day);
            end
        case 'qts'
            if isMetric
                yl = 'Temperature (Degrees Kelvin)';
            else
                yl = 'Temperature (Degrees Fahrenheit)';
                d = d*1.8 + 32;
            end
        case {'bhp', 'pressure'}
            if isMetric
                yl = 'Pressure (Pascal)';
            else
                yl = 'Pressure (Barsa)';
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