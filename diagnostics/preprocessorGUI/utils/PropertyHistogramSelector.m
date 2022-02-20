%%
% This class was taken from the PropertyDisplaySelector (mrst-visualization/mrst-gui/utils)
% and modified for the histogram selector in the EnsembleGUI. Some of the animation 
% funcionality of PropertyDisplaySelector are comented here and can be
% used in the future for dynamic histograms functionalities


classdef PropertyHistogramSelector < UIItem
    properties
        Callback
        propPopup
        minSlider = [];
        maxSlider = [];
        minEdit = [];
        maxEdit = [];
        binsSlider = [];
        binsEdit   = [];
        singleStep = false;
        renderTime   = .01;
        props
        enableSwitch = true;
        logSwitchBox
        statisticsForAll = false;
    end
    properties (Dependent)
        selection
        propIx
        statIx
        Min     % min prop value
        Max     % max prop value
        bins    % number of bins for the histograms
        minValue     % selected lower value
        maxValue     % selected higher value
        binsValue    % selected number of histobrams bins
        logSwitch
    end
    
    methods
        
        function s = PropertyHistogramSelector(varargin)
            % sampleprops must be updated with limits! will produce error
            sampleProps = struct('static',  struct('name', {{'sp1', 'sp2'}}, 'limits', ones(2,1)*[0 1]));
                             
            opt = struct('Parent',          [], ...
                         'Callback',        'disp(''Hello'')', ...
                         'Position', [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title', 'Select property for display', ...
                         'props', sampleProps, ...
                         'includeFilter', false, ...
                         'includeLogSwitch', false, ...
                         'statisticsForAll', false);
            [opt, extraOpt] = merge_options(opt, varargin{:});
             
            propPopup = uicontrol('Parent', [], 'Style', 'popup', 'String', opt.props.static.name, ...
                                  'Value', 1, 'Visible', 'off', 'UserData', opt.props);
            controls      = {{propPopup, []}};
            controlLayout = {[nan, .3]}; 
            
            if opt.includeLogSwitch
                logSwitchBox = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'log10', ...
                                         'Visible', 'off');
                controls{1}{2} =  logSwitchBox;
            end
                
            
            if opt.includeFilter
                minText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                      'String', 'Min ', 'Visible', 'off');
                maxText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                      'String', 'Max', 'Visible', 'off');
                binsText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                      'String', 'Bins', 'Visible', 'off');
                minSlider = uicontrol('Parent', [], 'Style', 'slider','Tag', 'min', ...
                                      'Min', 0, 'Max', 1, 'Value', 0);
                maxSlider = uicontrol('Parent', [], 'Style', 'slider','Tag', 'max',...
                                      'Min', 0, 'Max', 1, 'Value', 1); 
                binsSlider = uicontrol('Parent', [], 'Style', 'slider',...
                                      'Min', 1, 'Max', 500, 'Value', 80); 
                minEdit   = uicontrol('Parent', [], 'Style', 'edit', ... 
                                      'String', num2str(0), 'Min', 0, 'Max', 1);                            
                maxEdit   = uicontrol('Parent', [], 'Style', 'edit', ... 
                                      'String', num2str(1),'Min', 0, 'Max', 1);
                binsEdit   = uicontrol('Parent', [], 'Style', 'edit', ... 
                                      'String', num2str(80),'Min', 1, 'Max', 1);
                controls      = [controls, {{minText, minSlider, minEdit}, {maxText, maxSlider, maxEdit},{binsText,binsSlider,binsEdit}}];
                controlLayout = [controlLayout, {[.15, nan, .3], [.15, nan, .3],[.15, nan, .3]}];
            end
                        
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                         'controlWidths', controlLayout, 'Title', opt.Title, ...
                         'Position', opt.Position,'Visible', 'off', extraOpt{:});
            
            s.props = opt.props;                            
            s.propPopup = propPopup;
            s.fixedHeight = true;
            % main callback
            s.Callback = opt.Callback;
            % item callbacks (include main)
            propPopup.Callback = @s.propCallback;
                        
            if opt.includeLogSwitch
                s.logSwitchBox = logSwitchBox;
                s.logSwitchBox.Value  = false;
                logSwitchBox.Callback = @s.logSwitchCallback;
            end
           
            if opt.includeFilter
                s.minSlider  = minSlider;
                s.maxSlider  = maxSlider;
                s.binsSlider = binsSlider;
                s.minEdit    = minEdit;
                s.maxEdit    = maxEdit;
                s.binsEdit   = binsEdit;
                
                s.minSlider.Callback  = @s.sliderCallback;
                s.maxSlider.Callback  = @s.sliderCallback;
                s.binsSlider.Callback = @s.sliderCallback;
                s.minEdit.Callback    = @s.editCallback;
                s.maxEdit.Callback    = @s.editCallback;
                s.binsEdit.Callback   = @s.editCallback;
            end
            
        end
                
        function set.propIx(s, val)
            s.propPopup.Value = val;
        end
        function val = get.propIx(s)
            val = s.propPopup.Value;
        end
        
        function val = get.selection(s)
            val = s.propPopup.String{s.propIx};
        end
        
        function set.statIx(s, val)
            s.statPopup.Value = val;
        end
        function val = get.statIx(s)
            val = s.statPopup.Value;
        end  
        
        function set.Min(s, val)
            s.minValue = val;

            s.minSlider.Min = val;
            s.maxSlider.Min = val;
        end
        function val = get.Min(s)
            val = s.minSlider.Min; % should be the same as maxSlider
        end
        
        function set.Max(s, val)
            s.maxValue = val;

            s.minSlider.Max = val;
            s.maxSlider.Max = val;
        end
        function val = get.Max(s)
            val = s.maxSlider.Max; % should be the same as maxSlider
        end
        
        function set.bins(s, val)
            s.binsValue = round(val);
        end
        function val = get.bins(s)
            val = round(s.binsSlider.Value); % should be the same as binsSlider and it should be integer but I dont know where  should be specified as integer so I will round everything
        end
        
        function set.minValue(s, val)
            s.minEdit.String = num2str(val,'%.4g');

            s.minSlider.Value = val;
        end
        function val = get.minValue(s)
            val = s.minSlider.Value;
        end
        
        function set.maxValue(s, val)
            s.maxEdit.String = num2str(val,'%.4g');

            s.maxSlider.Value = val;
        end
        function val = get.maxValue(s)
            val = s.maxSlider.Value;
        end
        
        function set.binsValue(s, val)
            s.binsEdit.String = num2str(round(val),'%d');
            s.binsSlider.Value = round(val);
        end
        function val = get.binsValue(s)
            val = round(s.binsSlider.Value);
        end
        
        
        function set.logSwitch(s, val)
            s.logSwitchBox.Value = double(val);
        end
        function val = get.logSwitch(s)
            val = logical(s.logSwitchBox.Value);
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
            if ~isempty(s.props.static.limits)
                lims = s.props.static.limits{s.propIx};
                if s.logSwitch
                    lims = makeLogCompatible(lims);
                end
                s.Min = lims(1);
                s.Max = lims(2);
                s.bins= 80;
            end
        end
        
        function sliderCallback(s, src, event)
            switch src.Tag
                case 'max'
                    s.minValue  = s.minSlider.Value;
                    s.maxValue  = max(s.maxSlider.Value, s.minSlider.Value*1.01);
                    s.binsValue = s.binsSlider.Value;
                case 'min'
                    s.minValue  = min(s.minSlider.Value, s.maxSlider.Value*.99);
                    s.maxValue  = s.maxSlider.Value;
                    s.binsValue = s.binsSlider.Value;
                otherwise
                    s.minValue  = s.minSlider.Value;
                    s.maxValue  = s.maxSlider.Value;
                    s.binsValue = s.binsSlider.Value;
            end
            s.selectorCallback(src, event);
        end
        
        function editCallback(s, src, event)
            doCallback = false;
            if all(ismember(s.minEdit.String, '0123456789-+.eE'))
                val = max(s.Min, min(s.Max, str2double(s.minEdit.String)));
                s.minValue = min(val, s.maxValue*0.99);
                doCallback = true;
            end
            if all(ismember(s.maxEdit.String, '0123456789-+.eE'))
                val = max(s.Min, min(s.Max, str2double(s.maxEdit.String)));
                s.maxValue = max(val, s.minValue*1.01);
                doCallback = true;
            end
            if all(ismember(s.binsEdit.String, '0123456789-+.eE'))
                val = max(1, min(200, str2double(s.binsEdit.String)));
                s.binsValue = round(val); % To round just in case, TODO, check the recursivity of using round
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
        
        
        function logSwitchCallback(s, src, event)
            % change scale of slidebars if present, but keep selection
          %  if ~(isempty(s.minSlider)||isempty(s.maxSlider))
          %      s.updateSlideBars(src, event)
          %  end
            s.selectorCallback(src, event);
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
