classdef TrajectorySliceSelector < UIItem
    properties
        Callback
        angleSlider = [];
        angleEdit   = [];
        nPntPopup
        alphaEdit
        angleLimits = [-90, 90];
    end
    properties (Dependent)
        angle
        nPnt
        alpha
    end
    
    methods
        
        function s = TrajectorySliceSelector(varargin)                 
            opt = struct('Parent',          [], ...
                         'Callback',        @(src,event)disp('Hello'), ...
                         'Position', [10 10 300 300], ...
                         'Visible',         'on', ...
                         'Title', 'Trajectory slice options');
                     
            [opt, extraOpt] = merge_options(opt, varargin{:});
             
            angleText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                    'String', 'Slice angle ', 'Visible', 'off');
            angleSlider = uicontrol('Parent', [], 'Style', 'slider',...
                                    'Min', -90, 'Max', 90, 'Value', -90, 'Tag', 'angle', ...
                                    'SliderStep', [1/60 1/18]); 
            angleEdit   = uicontrol('Parent', [], 'Style', 'edit', ... 
                                    'String', num2str(-90), 'Min', 0, 'Max', 1, ...
                                    'Tag', 'angle'); 
                                
            nPntText    = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                    'String', '# points', 'Visible', 'off');
            sel = cellfun(@num2str, num2cell(2:10), 'UniformOutput', false);
            nPntPopup   = uicontrol('Parent', [], 'Style', 'popup', 'String', sel, ...
                                    'Value', 1, 'Visible', 'off', 'Tag', 'nPnt');
                                
            alphaText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                    'String', 'alpha ', 'Visible', 'off');
            alphaEdit   = uicontrol('Parent', [], 'Style', 'edit', ... 
                                    'String', '0.2', 'Min', 0, 'Max', 1, 'Tag', 'alpha'); 
            
                                
            controls      = {{angleText, angleSlider, angleEdit}, ...
                             {nPntText, nPntPopup, alphaText, alphaEdit}};
            controlLayout = {[.15, nan, nan], [.25, .25 .25 .25]}; 
            
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                         'controlWidths', controlLayout, 'Title', opt.Title, ...
                         'Position', opt.Position,'Visible', 'off', extraOpt{:});
            
            s.angleSlider = angleSlider;
            s.angleEdit   = angleEdit;
            s.nPntPopup   = nPntPopup;
            s.alphaEdit   = alphaEdit;
            s.fixedHeight = true;
            % main callback
            s.Callback = opt.Callback;
            % item callbacks (include main)
            %typePopup.Callback = @s.typeCallback;
            %propPopup.Callback = @s.propCallback;
            %statPopup.Callback = @s.statCallback;
            
            s.angleSlider.Callback = @s.angleSliderCallBack;
            s.angleEdit.Callback   = @s.angleEditCallBack;
            s.nPntPopup.Callback   = @s.angleSliderCallBack;
            s.alphaEdit.Callback   = @s.angleSliderCallBack;
           
            % set visible
            s.Visible = opt.Visible;
        end
        
        function set.angle(s, val)
            val = max(s.angleLimits(1), min(s.angleLimits(2), val));
            s.angleSlider.Value = val;
            s.angleEdit.String  = num2str(val);
        end
        function val = get.angle(s)
            val = s.angleSlider.Value;
        end
        
        function set.nPnt(s, val)
            ix = strcmp(num2str(val), s.nPntPopup.String);
            if nnz(ix) == 1
                s.nPntPopup.Value = find(ix);
            end
        end
        function val = get.nPnt(s)
            val = str2double(s.nPntPopup.String(s.nPntPopup.Value));
        end
        
        function set.alpha(s, val)
            val = max(0, min(1, val));
            s.alphaEdit.String = num2str(val);
        end
        function val = get.alpha(s)
            val = str2double(s.alphaEdit.String);
        end  
                
        function angleSliderCallBack(s, src, event)
            val = s.angleSlider.Value;
            s.angleEdit.String = num2str(val);
            s.Callback(src, event);
        end
        
        function angleEditCallBack(s, src, event)
            str = s.angleEdit.String;
            s.angle = str2double(str);
            s.Callback(src, event);
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
