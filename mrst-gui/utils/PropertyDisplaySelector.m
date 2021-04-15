classdef PropertyDisplaySelector < UIItem
    properties
        Callback
        typePopup
        propPopup
        statPopup
        minSlider = [];
        maxSlider = [];
        minEdit = [];
        maxEdit = [];
        singleStep = false;
        %startPlayFcn = [];
        %stopPlayFcn  = [];
        renderTime   = .01;
        playMode     = 'stop';
        playRange    = [];
        playDuration = 1;
        playBar      = [];
        props
        enableSwitch = true;
        logSwitchBox
        statisticsForAll = false;
    end
    properties (Dependent)
        typeIx
        propIx
        statIx
        Min     % min prop value
        Max     % max prop value
        minValue     % selected lower threshold
        maxValue    % selected higher threshold
        logSwitch
    end
    
    methods
        
        function s = PropertyDisplaySelector(varargin)
            % sampleprops must be updated with limits! will produce error
            sampleProps = struct('static',  struct('name', {{'sp1', 'sp2'}}, 'limits', ones(2,1)*[0 1]), ...
                                 'dynamic', struct('name', {{'dp1', 'dp2', 'dp3'}}, 'limits', ones(3,1)*[0 1]), ...
                                 'diagnostics', struct('name', {{'diagn1', 'diagn'}}, 'limits', ones(2,1)*[0 1]));
                             
            opt = struct('Parent',          [], ...
                         'Callback',        'disp(''Hello'')', ...
                         'Position', [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title', 'Select property for display', ...
                         'props', sampleProps, ...
                         'startPlayCallback', [], ...
                         'stopPlayCallback', [], ...
                         'includeFilter', false, ...
                         'includePlayer', false, ...
                         'includeEnableSwitch', false, ...
                         'includeLogSwitch', false, ...
                         'statisticsForAll', false);
            [opt, extraOpt] = merge_options(opt, varargin{:});
             
            typePopup = uicontrol('Parent', [], 'Style', 'popup', 'String', fieldnames(opt.props), ...
                                  'Value', 1, 'Visible', 'off');
            propPopup = uicontrol('Parent', [], 'Style', 'popup', 'String', opt.props.static.name, ...
                                  'Value', 1, 'Visible', 'off', 'UserData', opt.props);
            statPopup = uicontrol('Parent', [], 'Style', 'popup', 'String', {'n/a'}, ...
                                  'Value', 1, 'Visible', 'off', 'Enable', 'off');
            controls      = {{typePopup, []}, {propPopup, statPopup}};
            controlLayout = {[nan, .3],[nan, .3]}; 
            
            if opt.includeLogSwitch
                logSwitchBox = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'log10', ...
                                         'Visible', 'off');
                controls{1}{2} =  logSwitchBox;
            end
                
            if opt.includeEnableSwitch
                enableSwitch = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'Enable', ...
                                         'Visible', 'off');
                controls =  [{{enableSwitch}}, controls]; 
                controlLayout = [{nan}, controlLayout];
            end
            
            if opt.includeFilter
                minText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                      'String', 'Min ', 'Visible', 'off');
                maxText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                      'String', 'Max', 'Visible', 'off');
                minSlider = uicontrol('Parent', [], 'Style', 'slider',...
                                      'Min', 0, 'Max', 1, 'Value', 0);
                maxSlider = uicontrol('Parent', [], 'Style', 'slider',...
                                      'Min', 0, 'Max', 1, 'Value', 1);                             
                minEdit   = uicontrol('Parent', [], 'Style', 'edit', ... 
                                      'String', num2str(0), 'Min', 0, 'Max', 1);                            
                maxEdit   = uicontrol('Parent', [], 'Style', 'edit', ... 
                                      'String', num2str(1),'Min', 0, 'Max', 1);
                controls      = [controls, {{minText, minSlider, minEdit}, {maxText, maxSlider, maxEdit}}];
                controlLayout = [controlLayout, {[.15, nan, nan], [.15, nan, nan]}];
            end
            
            
            
            if opt.includePlayer
                if ~opt.includeFilter
                    warning('Player is not included unless optioin ''inludeFilter'' is set to true')
                end
                playReverseButton = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'playReverse', ...
                                              'CData', getIcon('ic_playreverse.png'));
                pauseButton       = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'pause', ...
                                              'CData', getIcon('ic_pause.png'));  
                playButton        = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'play', ...
                                              'CData', getIcon('ic_play.png'));    
                stopButton        = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'stop', ...
                                              'CData', getIcon('ic_stop.png')); 
                durationPopup     = uicontrol('Parent', [], 'Style', 'popup', ... 
                                              'String', {'1 sec', '2 sec', '5 sec', '10 sec', '20 sec'},...
                                              'Value', 1, 'UserData', [1 2 5 10 20]);
                controls      = [controls, {{playReverseButton, pauseButton, playButton, stopButton, durationPopup}}];
                controlLayout = [controlLayout, {[.15 .15 .15 .15 nan]}];
            end
            
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                         'controlWidths', controlLayout, 'Title', opt.Title, ...
                         'Position', opt.Position,'Visible', 'off', extraOpt{:});
            
            s.props = opt.props;                            
            s.typePopup = typePopup;
            s.propPopup = propPopup;
            s.statPopup = statPopup;
            s.fixedHeight = true;
            % main callback
            s.Callback = opt.Callback;
            % item callbacks (include main)
            typePopup.Callback = @s.typeCallback;
            propPopup.Callback = @s.propCallback;
            statPopup.Callback = @s.statCallback;
            
            if opt.statisticsForAll
                s.statisticsForAll = true;
            end
            
            if opt.includeEnableSwitch
                s.enableSwitch = false;
                enableSwitch.Callback = @s.eswitchCallback;
            end
            
            if opt.includeLogSwitch
                s.logSwitchBox = logSwitchBox;
                s.logSwitchBox.Value  = false;
                logSwitchBox.Callback = @s.logSwitchCallback;
            end
           
            if opt.includeFilter
                s.minSlider = minSlider;
                s.maxSlider = maxSlider;
                s.minEdit   = minEdit;
                s.maxEdit   = maxEdit;
                
                s.minSlider.Callback = @s.sliderCallback;
                s.maxSlider.Callback = @s.sliderCallback;
                s.minEdit.Callback   = @s.editCallback;
                s.maxEdit.Callback   = @s.editCallback;
            end
            
            if opt.includePlayer && opt.includeFilter
                s.playBar = [playReverseButton, pauseButton, playButton, stopButton, durationPopup];
                fnc = @(src, event)s.animButtonCallback(src, event, opt.startPlayCallback, opt.stopPlayCallback);
                playReverseButton.Callback = fnc;
                pauseButton.Callback       = fnc;
                playButton.Callback        = fnc;
                stopButton.Callback        = fnc;
                durationPopup.Callback     = @s.durationCallback;
            end
           
            % set visible
            s.Visible = opt.Visible;
            if opt.includeEnableSwitch
               s.Enable = 'off';
               enableSwitch.Enable = 'on';
            end
            
            s.setTypeEnable(s);
        end
        
        function set.typeIx(s, val)
            s.typePopup.Value = val;
            tp = s.typePopup.String{val};
            s.propPopup.String = s.propPopup.UserData.(tp).name;
        end
        function val = get.typeIx(s)
            val = s.typePopup.Value;
        end
        
        function set.propIx(s, val)
            s.propPopup.Value = val;
        end
        function val = get.propIx(s)
            val = s.propPopup.Value;
        end
        
        function set.statIx(s, val)
            s.statPopup.Value = val;
        end
        function val = get.statIx(s)
            val = s.statPopup.Value;
        end  
        
        function set.Min(s, val)
            s.minValue = val;
            if s.logSwitch
                val = log10(val);
            end
            s.minSlider.Min = val;
            s.maxSlider.Min = val;
        end
        function val = get.Min(s)
            val = s.minSlider.Min; % should be the same as maxSlider
            if s.logSwitch
                val = 10^val;
            end
        end
        
        function set.Max(s, val)
            s.maxValue = val;
            if s.logSwitch
                val = log10(val);
            end
            s.minSlider.Max = val;
            s.maxSlider.Max = val;
        end
        function val = get.Max(s)
            val = s.maxSlider.Max; % should be the same as maxSlider
            if s.logSwitch
               val = 10^val;
            end
        end
        
        function set.minValue(s, val)
            s.minEdit.String = num2str(val,'%7.2g');
            if s.logSwitch
               val = log10(val);
            end
            s.minSlider.Value = val;
        end
        function val = get.minValue(s)
            val = s.minSlider.Value;
            if s.logSwitch
                val = 10^val;
            end
        end
        
        function set.maxValue(s, val)
            s.maxEdit.String = num2str(val,'%7.2g');
            if s.logSwitch
                val = log10(val);
            end
            s.maxSlider.Value = val;
        end
        function val = get.maxValue(s)
            val = s.maxSlider.Value;
            if s.logSwitch
                val = 10^val;
            end
        end
        
        function set.logSwitch(s, val)
            s.logSwitchBox.Value = double(val);
        end
        function val = get.logSwitch(s)
            val = logical(s.logSwitchBox.Value);
        end
            
        function typeCallback(s, src, event)
            s.typeIx = s.typePopup.Value;
            s.setTypeEnable(s);
           tp = s.typePopup.String{s.typeIx};
           s.propPopup.Value  = 1;
           s.propPopup.String = s.propPopup.UserData.(tp).name;
           if ~isempty(s.minSlider)
               s.updateSlideBars(src, event);
           end
           s.selectorCallback(src, event)
        end
        
        function propCallback(s, src, event)
            % if s.propIx ~= s.propPopup.Value
            s.propIx = s.propPopup.Value;
            if ~isempty(s.minSlider)
                s.updateSlideBars(src, event);
            end
            s.selectorCallback(src, event)
        end
        
        function statCallback(s, src, event)
            s.statIx = s.statPopup.Value;
            s.selectorCallback(src, event);
        end
        
        function updateSlideBars(s, src, event)
            tp = s.typePopup.String{s.typeIx};
            if ~isempty(s.props.(tp).limits)
                lims = s.props.(tp).limits{s.propIx};
                if s.logSwitch
                    lims = makeLogCompatible(lims);
                end
                s.Min = lims(1);
                s.Max = lims(2);
            end
        end
        
        function sliderCallback(s, src, event)
            if ~s.logSwitch
                s.minValue = s.minSlider.Value;
                s.maxValue = s.maxSlider.Value;
            else
                s.minValue = 10^s.minSlider.Value;
                s.maxValue = 10^s.maxSlider.Value;
            end
            s.selectorCallback(src, event);
        end
        
        function editCallback(s, src, event)
            doCallback = false;
            if all(ismember(s.minEdit.String, '0123456789-+.eE'))
                val = max(s.Min, min(s.Max, str2double(s.minEdit.String)));
                s.minValue = val;
                doCallback = true;
            end
            if all(ismember(s.maxEdit.String, '0123456789-+.eE'))
                val = max(s.Min, min(s.Max, str2double(s.maxEdit.String)));
                s.maxValue = val;
                doCallback = true;
            end
            if doCallback
                s.selectorCallback(src, event);
            end
        end
        
        function selectorCallback(s, src, event)
            %main callback
            if ischar(s.Callback)
                eval(s.Callback);
            else
                s.Callback(src, event);
            end
        end
        function eswitchCallback(s, src, event)
            if src.Value == 1
                s.Enable = 'on';
            else
                s.Enable = 'off';
                src.Enable = 'on';
            end
            s.selectorCallback(src, event);
        end
        
        function logSwitchCallback(s, src, event)
            % change scale of slidebars if present, but keep selection
            if ~(isempty(s.minSlider)||isempty(s.maxSlider))
                s.updateSlideBars(src, event)
            end
            s.selectorCallback(src, event);
        end
        
        function durationCallback(s, src, ~)
            s.playDuration = src.UserData(src.Value);
        end
        
        function animButtonCallback(s, src, event, startPlay, stopPlay)
            instr = src.Tag;
            if ~strcmp(instr, s.playMode)   % instruction different from current mode
                if strcmp(instr, 'stop')    % stop playmode
                    s.playMode = instr;
                    if ~isempty(stopPlay)   % run stop callback
                        stopPlay(src, event);
                    end
                    [s.minValue, s.maxValue] = deal(s.playRange(1), s.playRange(2));
                    s.selectorCallback(src, event);
                else
                    if ~( strcmp(s.playMode, 'stop') && strcmp(instr, 'pause') )
                        curMax = s.maxValue;
                        if strcmp(s.playMode, 'stop')   % start playmode
                            if ~isempty(startPlay)
                                startPlay(src, event);  % run start callback
                            end
                            s.playRange = [s.minValue, s.maxValue];
                            if strcmp(instr, 'play')
                                curMax = s.minValue;
                            end
                        end
                        s.playMode  = instr;
                        
                        while ~any(strcmp(s.playMode, {'stop', 'pause'})) && strcmp(instr, s.playMode)
                            [inc, ps] = getPlayIncrement(s);
                            if ~s.logSwitch
                                curMax = curMax + inc;
                                % rewind
                                curMax = mod(curMax - s.playRange(1), ...
                                         diff(s.playRange)) + s.playRange(1);
                            else
                                logCurMax = log10(curMax)+inc;
                                % rewind
                                logCurMax = mod(logCurMax - log10(s.playRange(1)), ...
                                         diff(log10(s.playRange))) + log10(s.playRange(1));
                                curMax = 10^(logCurMax);
                            end
                            % rewind
                            curMax = mod(curMax - s.playRange(1), ...
                                diff(s.playRange)) + s.playRange(1);
                            s.maxValue  = curMax;
                            st = tic;
                            s.Callback(src, event);
                            drawnow limitrate;
                            pause(ps);
                            s.renderTime = .8*s.renderTime + .2*toc(st);
                        end
                    end
                end
            end
        end 
            
    end
    methods (Static)
        function setTypeEnable(s)
            curType  = s.typePopup.String{s.typeIx};
            if (strcmp(curType, 'static') && ~s.statisticsForAll) || s.singleStep
                s.statPopup.Enable = 'off';
                s.statPopup.Value  = 1;
                s.statPopup.String = 'n/a';
            else
                s.statPopup.Enable = 'on';
                s.statPopup.String = {'mean', 'std', 'max diff'};
            end
        end
    end
end

% -------------------------------------------------------------------------

function lims = makeLogCompatible(lims)
if lims(2) <=0
    error('All data negative, can''t take log')
else
    lims(1) = max(lims(1), lims(2)*10^(-5));
end
end

% -------------------------------------------------------------------------

function [inc, ps] = getPlayIncrement(s)
dur     = s.playDuration;
nframes = min(dur*50, dur/s.renderTime);
sgn = 1;
if strcmp(s.playMode, 'playReverse')
    sgn = -1;
end
if ~s.logSwitch
    inc = sgn*diff(s.playRange)/nframes;
else
    inc = sgn*diff(log10(s.playRange))/nframes;
end
ps  = dur/nframes - s.renderTime;
ps = max(ps, 1/20);
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
