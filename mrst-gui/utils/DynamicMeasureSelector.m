classdef DynamicMeasureSelector < UIItem
    properties
        Callback
        leftPopup
        rightPopup
        leftAvgSwitch
        rightAvgSwitch
        props
        enableSwitch = true;
        panelNo
    end
    properties (Dependent)
        leftIx
        rightIx
        leftSwitch
        rightSwitch
    end
    
    methods
        
        function s = DynamicMeasureSelector(varargin)
            sampleProps = {{{'none',  'measure1', 'measure2'}}};
            opt = struct('Parent',          [], ...
                         'Callback',        'disp(''Hello'')', ...
                         'Position',        [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title',           'Heterogeneity measures', ...
                         'props',           sampleProps, ...
                         'includeEnableSwitch', false, ...
                         'includeAvgSwitch', false);
            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            leftText  = uicontrol('Parent', [], 'Style', 'text', ...
                                  'Value',  [], 'String', 'Left', 'Visible', 'off');
            leftPopup = uicontrol('Parent', [], 'Style', 'popup', ...
                                  'Value',  1, 'String', opt.props{:}, 'Visible', 'off');
                               
            rightText  = uicontrol('Parent', [], 'Style', 'text', ...
                                  'Value',  [], 'String', 'Right', 'Visible', 'off');
            rightPopup = uicontrol('Parent', [], 'Style', 'popup', ...
                                  'Value',  1, 'String', opt.props{:}, 'Visible', 'off');
                               
            controls      = {{leftText, leftPopup}, {rightText, rightPopup}};
            controlLayout = {[.2, nan],[.2, nan]}; 
            
            if opt.includeAvgSwitch
                leftAvgSwitch  = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'Avg', ...
                    'Visible', 'off');
                rightAvgSwitch = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'Avg', ...
                    'Visible', 'off');
                controls{1} = [controls{1}, {leftAvgSwitch}];
                controls{2} = [controls{2}, {rightAvgSwitch}];
                controlLayout = {[.2, nan, .25],[.2, nan, .25]}; 
            end
            
            if opt.includeEnableSwitch
                enableSwitch = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'Enable', ...
                    'Visible', 'off');
                controls =  [{{enableSwitch}}, controls]; 
                controlLayout = [{nan}, controlLayout];
            end
            
               
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                                'controlWidths', controlLayout, 'Title', opt.Title, ...
                                'Position', opt.Position,'Visible', 'off', extraOpt{:});
            
            s.props = opt.props;                            
            s.leftPopup      = leftPopup;
            s.rightPopup     = rightPopup;
            s.fixedHeight = true;
            % main callback
            s.Callback = opt.Callback;
            % item callbacks (include main)
            leftPopup.Callback      = @s.leftCallback;
            rightPopup.Callback     = @s.rightCallback;
            
            
            if opt.includeAvgSwitch
                s.leftAvgSwitch  = leftAvgSwitch;
                s.rightAvgSwitch = rightAvgSwitch;
                leftAvgSwitch.Callback  = @s.leftCallback;
                rightAvgSwitch.Callback = @s.rightCallback;
            end
            
            
            if opt.includeEnableSwitch
                s.enableSwitch = false;
                enableSwitch.Callback = @s.eswitchCallback;
            end
            
            % set visible
            s.Visible = opt.Visible;
            s.Enable = 'off';
            if opt.includeEnableSwitch
               s.Enable = 'off';
               enableSwitch.Enable = 'on';
            end
            
        end
        
        function set.leftIx(s, val)
            s.leftPopup.Value = val;
        end
        function val = get.leftIx(s)
            val = s.leftPopup.Value;
        end
        
        function set.rightIx(s, val)
            s.rightPopup.Value = val;
        end
        function val = get.rightIx(s)
            val = s.rightPopup.Value;
        end
        
        function set.leftSwitch(s, val)
            s.leftAvgSwitch.Value = double(val);
        end
        function val = get.leftSwitch(s)
            val = false;
            if ~isempty(s.leftAvgSwitch)
                val = logical(s.leftAvgSwitch.Value);
            end
        end
        
        function set.rightSwitch(s, val)
            s.rightAvgSwitch.Value = double(val);
        end
        function val = get.rightSwitch(s)
            val = false;
            if ~isempty(s.rightAvgSwitch)
                val = logical(s.rightAvgSwitch.Value);
            end
        end
        
        function leftCallback(s, src, event)
            %s.leftIx = s.leftPopup.Value;
            s.panelNo = 1;
            s.selectorCallback(src, event);
        end
        
        function rightCallback(s, src, event)
            %s.rightIx = s.rightPopup.Value;
            s.panelNo = 2;
            s.selectorCallback(src, event);
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
    end
    methods (Static)
        function setTypeEnable(s)
            curType  = s.typePopup.String{s.typeIx};
            if strcmp(curType, 'static') || s.singleStep
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
