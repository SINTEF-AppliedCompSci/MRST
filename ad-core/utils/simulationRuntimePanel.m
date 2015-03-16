function ok = simulationRuntimePanel(model, ctrl_reports, schedule, simtime, varargin)
    persistent simAbort simPause simDebug simPlaySound
    if isempty(simAbort)
        simAbort = false;
        simPause = false;
    end
    pause(0.01);
    
    opt = struct('figure', []);
    opt = merge_options(opt, varargin{:});
    
    
    ctrl_reports = ctrl_reports(cellfun(@(x) ~isempty(x), ctrl_reports));
    stepNo = numel(ctrl_reports);
    
    panelExists = ~isempty(findobj('Tag', 'mrst-simpanel'));
    if ~panelExists
        [simAbort, simPause, simPlaySound] = deal(false);
    end
    ok = true;

    % Redraw plot
    p = getPanel();
    
    
    txt = ['Simulating ', class(model), ' '];
    uicontrol(p, 'units', 'Normalized', ...
    'Style', 'text', 'string', txt,...
    'fontsize', 18, ...
    'position', [0.1, 0.9, .8, .1])

    
    % Control step stuff
    
    stepTotal = numel(schedule.step.val);
    done = stepNo/stepTotal;
    
    completed = stepNo == stepTotal;
    if completed
        str = ['Simulation done: ', num2str(stepTotal), ' control steps solved'];
    else
        str = ['Solving control step number ', num2str(stepNo + 1), ' of ', num2str(stepTotal)];
    end
    simpleUIbar(p, done, 0.75, .15, str)
    
    % Actual time
    t = sum(schedule.step.val(1:stepNo));
    T = sum(schedule.step.val);
    tm = @(x) lower(formatTimeRange(x));
    str = ['Simulated ', tm(t), ' of ', tm(T), ' in schedule'];
    simpleUIbar(p, t/T, 0.60, .15, str)
    
    % Simtime business
    simT = max(sum(simtime), 1);
    endT = max((stepTotal - stepNo)*simT/stepNo, 1);
    
    txt = ['Simulation has been running for ', tm(simT), ' '];
    if completed
        txt = [txt, 'and is done'];
    else
        txt = [txt, 'and will be completed in about ' tm(endT)];
    end
    txt = {txt, getStats(model)};
    
    uicontrol(p, 'units', 'Normalized', ...
        'Style', 'text', 'string', txt,...
        'position', [0.1, 0.5, .8, .1])

    axes('position', [0.1, 0.25, .8, .25])
    tmp = nan(stepTotal, 1);
    its = cellfun(@(x) x.Iterations, ctrl_reports);
    tmp(1:stepNo) = its;
    stairs(1:stepTotal, tmp);
    ylim([1, max(tmp) + 0.1])
    xlim([1, stepTotal])
    grid on
    
    txt = sprintf('Total number of iterations %d with an average of %1.2f iterations per control step', sum(its), sum(its)/stepNo);
    
    uicontrol(p, 'units', 'Normalized', ...
        'Style', 'text', 'string', txt,...
        'position', [0.1, 0.10, .8, .05])
    
    function pauseHandler(src, event)
        if get(src, 'Value');
            pause(.1)
            pauseHandler(src, event);
        end
    end

    function dbstopHandler(src, event)
        simDebug = get(src, 'Value');
    end

    function soundHandler(src, event)
        simPlaySound = get(src, 'Value');
    end
    
    function abortHandler(src, event)
        simAbort = true;
    end
    
    p2 = uipanel(p, 'position', [0 0 1 .10]);
    uicontrol(p2, 'units', 'Normalized', ...
        'Tag', 'cancelsim', ...
        'callback', @abortHandler, ...
        'Style', 'pushbutton', 'string', 'Abort',...
        'position', [0.05, .05, .3, .4])
    
    uicontrol(p2, 'units', 'Normalized', ...
        'Tag', 'playsound', 'value', simPlaySound,...
        'callback', @soundHandler, ...
        'Style', 'togglebutton', 'string', 'Play sound when done',...
        'position', [0.35, .05, .3, .4])
    
    uicontrol(p2, 'units', 'Normalized', ...
        'Tag', 'playsound', 'value', false,...
        'callback', @pauseHandler, ...
        'Style', 'togglebutton', 'string', 'Pause',...
        'position', [0.65, .05, .3, .4])
    
    if simPlaySound && completed
        d = load('gong.mat');
        soundsc(d.y);
    end
    
    
    if simAbort
        simAbort = false;
        ok = false;
        return
    end
    
    if simDebug
        dbstop
    end
    
    
    drawnow
end

function p = createPanel()
    try
        c = get(0,'defaultUicontrolBackgroundColor');
    catch
        c = [1 1 1];
    end
    df = get(0, 'defaultFigurePosition');
    p = figure('Tag', 'mrst-simpanel', 'Color', c, 'position', df.*[1 1 1 1.5]);
end

function h = getPanel()
    h = findobj(0, 'Tag', 'mrst-simpanel');
    if isempty(h)
        h = createPanel();
    end
    clf(h);
end

function simpleUIbar(parent, data, start, height, txt)
    h = height/2;
    set(0, 'CurrentFigure', parent);
    axes('position', [0.1, start + h, .8, h])
    plotProgress(data)
    
    uicontrol(parent, 'units', 'Normalized', ...
        'Style', 'text', 'string', txt,...
        'position', [0.1, start, .8, h*0.9])
end

function plotProgress(data)
    rectangle('Position', [0, 0.01, data, 0.99], 'FaceColor', 'r')
    hold on
    rectangle('Position', [0, 0.01, 1, 0.99], 'FaceColor', 'none')
    ylim([0, 1]);
    xlim([0, 1]);
    axis off
    txt = sprintf('%2.1f%%', 100*data);
    text(0.5, 0.5, txt)
end

function txt = getStats(model)
    txt = '';
    G = model.G;
    if ~isempty(G)
        txt = [txt, 'Grid contains ', num2str(G.cells.num), ' active cells. '];
        if ~isempty(model.rock)
            nr = size(model.rock.perm, 2);
            if nr == 1
                type = 'scalar';
            elseif (nr == 2 && G.griddim == 2) || (nr == 3 && G.griddim == 3)
                type = 'diagonal tensor';
            elseif (nr == 3 && G.griddim == 2) || (nr == 6 && G.griddim == 3)
                type = 'full tensor';
            else
                type = 'unknown';
            end
            txt = [txt, 'Model uses ', type, ' permeability. '];
        end
    end
end


