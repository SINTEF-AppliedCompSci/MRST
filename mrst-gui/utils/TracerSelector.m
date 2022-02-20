classdef TracerSelector < UIItem
    properties
        Callback
        leftPopup
        rightPopup
        watCheckBox
        totCheckBox
        regCheckBox
        nameEdit
        timeEdit
        props
        enableSwitch = true;
        panelNo
    end
    properties (Dependent)
        leftIx
        rightIx
        xAxisIx
        extendTime
    end
    
    methods
        function s = TracerSelector(varargin)
            sampleProps = {{'none', 'Estimate from TOF', 'Simulated'}};
            opt = struct('Parent',          [], ...
                         'Callback',        'disp(''Hello'')', ...
                         'Position',        [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title',           'RTD Distributions', ...
                         'props',           sampleProps, ...
                         'includeEnableSwitch', false);
            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            leftText  = uicontrol('Parent', [], 'Style', 'text', ...
                                  'Value',  [], 'String', 'Left', 'Visible', 'off');
            leftPopup = uicontrol('Parent', [], 'Style', 'popup', ...
                                  'Value',  1, 'String', opt.props, 'Visible', 'off');
                               
            rightText  = uicontrol('Parent', [], 'Style', 'text', ...
                                  'Value',  [], 'String', 'Right', 'Visible', 'off');
            rightPopup = uicontrol('Parent', [], 'Style', 'popup', ...
                                  'Value',  1, 'String', opt.props, 'Visible', 'off');
                              
            %xAxisText  = uicontrol('Parent', [], 'Style', 'text', ...
            %                      'Value',  [], 'String', 'x-axes', 'Visible', 'off');
            %xAxisPopup = uicontrol('Parent', [], 'Style', 'popup', ...
            %                       'Value',  1, 'String', {'time', 'cum production'}, 'Visible', 'off'); 
%             distText   = uicontrol('Parent', [], 'Style', 'text', ...
%                                    'Value',  [], 'String', 'Distributions -----', 'Visible', 'off');
%             watCheckBox = uicontrol('Parent', [], 'Style', 'checkbox', ...
%                                     'String', 'Water', 'Visible', 'off');
%             totCheckBox = uicontrol('Parent', [], 'Style', 'checkbox', ...
%                                     'String', 'Total', 'Visible', 'off');                        
            
%             regText     = uicontrol('Parent', [], 'Style', 'text', ...
%                                     'Value',  [], 'String', 'Cummulative 3D -----', 'Visible', 'off');
%             regCheckBox  = uicontrol('Parent', [], 'Style', 'checkbox', ...
%                                         'String', 'Include', 'Visible', 'off');
%             nameEdit     = uicontrol('Parent', [], 'Style', 'edit', ... 
%                                          'String', 'prop1', 'Min', 0, 'Max', 1);  
            
                               
            timeText = uicontrol('Parent', [], 'Style', 'text', ...
                                   'Value',  [], 'String', 'Time span (years):', 'Visible', 'off');
            timeEdit = uicontrol('Parent', [], 'Style', 'edit', ... 
                                   'String', num2str(20), 'Min', 0, 'Max', 1);  
                               
            %computeButton = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'compute', ...
            %                          'String', 'Compute');
                               
            controls      = {{leftText, leftPopup}, {rightText, rightPopup}, ...
                             {timeText, timeEdit}};
            controlLayout = {[.2, nan],[.2, nan], [.8 .2]}; 
            
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
            s.leftPopup = leftPopup;
            s.rightPopup = rightPopup;
            %s.xAxisPopup = xAxisPopup;
%            s.watCheckBox = watCheckBox;
 %           s.totCheckBox = totCheckBox;
%            s.regCheckBox = regCheckBox;
            s.timeEdit = timeEdit;
            
            s.timeEdit = timeEdit;
            s.fixedHeight = true;
            % main callback
            s.Callback = opt.Callback;
            % item callbacks (include main)
            leftPopup.Callback = @s.leftCallback;
            rightPopup.Callback = @s.rightCallback;
            %xAxisPopup.Callback = @s.xAxisCallback;
            timeEdit.Callback = @s.extendCallback;
            %computeButton.Callback = @s.Callback;
            
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
        
        function set.xAxisIx(s, val)
            s.xAxisPopup.Value = val;
        end
        function val = get.xAxisIx(s)
            val = s.xAxisPopup.Value;
        end
        
        function set.extendTime(s, val)
            s.timeEdit.String = num2str(val);
        end
        function val = get.extendTime(s)
            val = str2double(s.timeEdit.String);
        end
        
        function leftCallback(s, src, event)
            s.leftIx = s.leftPopup.Value;
            s.panelNo = 1;
            s.selectorCallback(src, event);
        end
        
        function rightCallback(s, src, event)
            s.rightIx = s.rightPopup.Value;
            s.panelNo = 2;
            s.selectorCallback(src, event);
        end
                
        function xAxisCallback(s, src, event)
            s.xAxisIx = s.xAxisPopup.Value;
            s.selectorCallback(src, event);
        end
        
        function extendCallback(s, src, event)
            s.extendTime = str2double(s.timeEdit.String);
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
