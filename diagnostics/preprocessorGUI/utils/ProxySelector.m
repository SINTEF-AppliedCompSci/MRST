classdef ProxySelector < UIItem
    properties
        Callback
        objectiveButton
        gradientButton
        slider
        controlButton
    end
    
    methods
        
        function s = ProxySelector(varargin)                 
            opt = struct('Parent',          [], ...
                         'Callback',        @(src,event)disp('Hello'), ...
                         'Position', [10 10 300 300], ...
                         'Visible',         'on', ...
                         'Title', 'Proxy and trajectory gradients');
                     
            [opt, extraOpt] = merge_options(opt, varargin{:});
            getui = @(style, string, value, tag)uicontrol('Parent', [], 'Style', style, 'Value', value, ...
                    'Visible', 'off', 'String', string, 'Tag', tag);
            getuitext = @(string)getui('text', string, [], '');
            
             
            objText      = getuitext('Case objectives:');
            objButton    = getui('pushbutton', 'Compute/plot', [], 'objective');
            gradText     = getuitext('Trajectory gradient:');
            gradButton   = getui('pushbutton', 'Compute/show', [], 'gradient');
            sliderText   = getuitext('Move along gradient:');
            slider       = uicontrol('Parent', [], 'Style', 'slider',...
                                    'Min', 0, 'Max', 1, 'Value', 0, 'Tag', 'slider', ...
                                    'SliderStep', [1/20 1/10]); 
            controlText     = getuitext('Control gradient:');
            controlButton   = getui('pushbutton', 'Compute/show', [], 'control');
                                 
            controls      = {{objText, objButton}, ...
                             {gradText, gradButton}, ...
                             {sliderText, slider}, ...
                             {controlText, controlButton}};
            controlLayout = {[.5 .5], [.5 .5], [.5 .5], [.5 .5]}; 
            
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                         'controlWidths', controlLayout, 'Title', opt.Title, ...
                         'Position', opt.Position,'Visible', 'off', extraOpt{:});
            
            s.objectiveButton  = objButton;
            s.gradientButton   = gradButton;
            s.slider           = slider;
            s.controlButton     = controlButton;
            s.fixedHeight = true;
            % main callback
            s.Callback = opt.Callback;
            s.objectiveButton.Callback = @s.Callback;
            s.gradientButton.Callback  = @s.gradientCallback;
            s.slider.Callback          = @s.Callback;
            s.controlButton.Callback   = @s.Callback;
            s.slider.Enable = 'off';
            % set visible
            s.Visible = opt.Visible;
        end
        
        function gradientCallback(s, src, event)
            s.Callback(src, event);
            s.slider.Enable = 'on';
            s.slider.Value  = 0;
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
