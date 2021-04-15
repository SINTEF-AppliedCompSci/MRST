classdef LineplotSelector < UIItem
    properties
        Callback
        popup
        props
        xLim
    end
    properties (Dependent)
        ix
        selection
        xLims
        yLims
    end
    
    methods
        
        function s = LineplotSelector(varargin)
            
            sampleProps = {{{'none',  'v2', 'v2'}}};
            opt = struct('Parent',          [], ...
                         'Callback',        @(src, event)disp('Boo'), ...
                         'Position',        [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title',           'Line plot', ...
                         'props',           sampleProps, ...
                         'includeXLims',    false, ...
                         'includeYLims',    false);
            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            %xText  = uicontrol('Parent', [], 'Style', 'text', ...
            %                      'Value',  [], 'String', 'x-axis', 'Visible', 'off');
            popup = uicontrol('Parent', [], 'Style', 'popup', ...
                              'Value',  1, 'String', opt.props{:}, 'Visible', 'off');
                               
            yText  = uicontrol('Parent', [], 'Style', 'text', ...
                                  'Value',  [], 'String', 'y-axis', 'Visible', 'off');
            yPopup = uicontrol('Parent', [], 'Style', 'popup', ...
                                  'Value',  1, 'String', opt.props{:}, 'Visible', 'off');
                               
            controls      = {{xText, xPopup}, {yText, yPopup}};
            controlLayout = {[.2, nan],[.2, nan]}; 
            
            if opt.includeFreezeSwitch
                xFreezeSwitch  = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'freeze', ...
                    'Tag', 'xFreeze', 'Visible', 'off');
                yFreezeSwitch = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'freeze', ...
                    'Tag', 'yFreeze', 'Visible', 'off');
                controls{1} = [controls{1}, {xFreezeSwitch}];
                controls{2} = [controls{2}, {yFreezeSwitch}];
                controlLayout = {[.2, nan, .25],[.2, nan, .25]}; 
            end
               
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                                'controlWidths', controlLayout, 'Title', opt.Title, ...
                                'Position', opt.Position,'Visible', 'off', extraOpt{:});
            
            s.props = opt.props;                            
            s.xPopup      = xPopup;
            s.yPopup     = yPopup;
            s.fixedHeight = true;
            % main callback
            s.Callback = opt.Callback;
            % item callbacks (include main)
            xPopup.Callback      = @s.Callback;
            yPopup.Callback      = @s.Callback;
            
            
            if isempty(opt.axes) || ~isvalid(opt.axes)
                s.Axes = [];
            else
                s.Axes = opt.axes;
            end
            
            if opt.includeFreezeSwitch
                s.xFreezeSwitch = xFreezeSwitch;
                s.yFreezeSwitch = yFreezeSwitch;
                xFreezeSwitch.Callback  = @s.Callback;
                yFreezeSwitch.Callback =  @s.Callback;
            end
            
            
            % set visible
            s.Visible = opt.Visible;
        end
        
        function set.xIx(s, val)
            s.xPopup.Value = val;
        end
        function val = get.xIx(s)
            val = s.xPopup.Value;
        end
        
        function val = get.xSelection(s)
            val = s.xPopup.String{s.xIx};
        end
        
        function set.yIx(s, val)
            s.yPopup.Value = val;
        end
        function val = get.yIx(s)
            val = s.yPopup.Value;
        end
        
        function val = get.ySelection(s)
            val = s.yPopup.String{s.yIx};
        end
        
        function set.xFreeze(s, val)
            s.xFreezeSwitch.Value = double(val);
            s.setAxesScale(s.xFreezeSwitch)
        end
        function val = get.xFreeze(s)
            val = logical(s.xFreezeSwitch.Value);
        end
        
        function set.yFreeze(s, val)
            s.yFreezeSwitch.Value = double(val);
            s.setAxesScale(s.yFreezeSwitch);
        end
        function val = get.yFreeze(s)
            val = logical(s.yFreezeSwitch.Value);
        end
   
%         function setAxesScale(s, src, ~)
%             if ~isempty(s.Axes) && isvalid(s.Axes)
%                 if src.Value == 1
%                     val = 'log';
%                 else
%                     val = 'linear';
%                 end
%                 if src.Tag == 'X'
%                     s.Axes.XScale = val;
%                 else
%                     s.Axes.YScale = val;
%                 end
%             end
%         end
            
%         function leftCallback(s, src, event)
%             %s.leftIx = s.leftPopup.Value;
%             s.panelNo = 1;
%             s.selectorCallback(src, event);
%         end
%         
%         function rightCallback(s, src, event)
%             %s.xIx = s.xPopup.Value;
%             s.panelNo = 2;
%             s.selectorCallback(src, event);
%         end
        
        function selectorCallback(s, src, event)
            %main callback
            if ischar(s.Callback)
                eval(s.Callback);
            else
                s.Callback(src, event);
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
