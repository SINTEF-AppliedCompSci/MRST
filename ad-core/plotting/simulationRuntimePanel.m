function ok = simulationRuntimePanel(model, states, ctrl_reports, solver, schedule, simtime, varargin)
% Internal function for drawing panel during simulation. See
% getPlotAfterStep.

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

    persistent simAbort simDebug simPlaySound simDump
    pause(0.01);
    
    opt = struct('figure', []);
    opt = merge_options(opt, varargin{:});
    
    
    ctrl_reports = ctrl_reports(cellfun(@(x) ~isempty(x), ctrl_reports));
    stepNo = numel(ctrl_reports);
    
    panelExists = ~isempty(findobj('Tag', 'mrst-simpanel'));
    if isempty(simAbort) || ~panelExists
        [simAbort, simDebug, simPlaySound, simDump] = deal(false);
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
    ylim([1, max(max(tmp) + 0.1, 1.01)])
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

    function dumpHandler(src, event)
        simDump = true;
    end
    
    function timestepHandler(src, event, factor)
        dt0 = solver.timeStepSelector.maxTimestep;
        if isinf(dt0)
            dt_new = 30*day;
        else
            dt_new = dt0*factor;
        end
        
        fprintf('**** Adjusted timestep interactively (dt: %s -> %s)\n', ...
                            formatTimeRange(dt0), formatTimeRange(dt_new));
        
        solver.timeStepSelector.maxTimestep = dt_new;
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
        'value', false,...
        'callback', @pauseHandler, ...
        'Style', 'togglebutton', 'string', 'Pause',...
        'position', [0.65, .05, .3, .4])
    
    uicontrol(p2, 'units', 'Normalized', ...
        'Tag', 'dbstop', ...
        'callback', @dbstopHandler, ...
        'Style', 'pushbutton', 'string', 'Debug',...
        'position', [0.05, .55, .3, .4])
    
    uicontrol(p2, 'units', 'Normalized', ...
        'callback', @dumpHandler, ...
        'Style', 'pushbutton', 'string', 'Dump to workspace',...
        'position', [0.35, .55, .3, .4])
    
    uicontrol(p2, 'units', 'Normalized', ...
        'callback', @(src, event) timestepHandler(src, event, 2), ...
        'Style', 'pushbutton', 'string', 'dtMax*2',...
        'position', [0.65, .55, .15, .4])
    
    uicontrol(p2, 'units', 'Normalized', ...
        'callback', @(src, event) timestepHandler(src, event, 1/2), ...
        'Style', 'pushbutton', 'string', 'dtMax/2',...
        'position', [0.8, .55, .15, .4])
    
    if simPlaySound && completed
        d = load('gong.mat');
        soundsc(d.y);
    end
    
    
    if simAbort
        simAbort = false;
        ok = false;
        return
    end
    
    if simDump
        data = struct();
        data.states = states;
        data.simtime = simtime;
        data.reports = ctrl_reports;
        assignin('base', 'simdata', data);
        simDump = false;
    end
    
    if simDebug
       stck = dbstack();
       c = stck(1);
       ln = sprintf('%d', c(1).line + 6);
       dbstop('in', c(1).file, 'at', ln)
       simDebug = false;
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
    pos = df.*[1 1 1 1.5];
    ssz = get(0,'ScreenSize');
    pos(3:4) = min([pos(3:4); ssz(3:4)-[0 90]]);
    pos(1:2) = min([pos(1:2); ssz(3:4)-pos(3:4)-[0 90]]);
    p = figure('Tag', 'mrst-simpanel', 'Color', c, 'position', pos);
end

function h = getPanel()
    h = findobj(0, 'Tag', 'mrst-simpanel');
    if isempty(h)
        h = createPanel();
    end
    clf(h);
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


