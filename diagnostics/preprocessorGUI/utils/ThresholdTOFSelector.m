classdef ThresholdTOFSelector < UIItem
    properties
        Callback
        tofSlider
        yearEdit
        pviEdit
        yearButton
        pviButton
        yearLimits
        pviLimits
    end
    properties (Dependent)
        Value
        timeUnit
    end
    
    methods
        
        function s = ThresholdTOFSelector(varargin)
            opt = struct('Parent',          [], ...
                         'Callback',        @(src, event)disp('Boo'), ...
                         'Position',        [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title',           'Threshold maximal time-of-flight', ...
                         'yearLimits',      [], ...
                         'pviLimits',       []);

            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            if isempty(opt.yearLimits), opt.yearLimits = [1 1000]; end
            if isempty(opt.pviLimits), opt.pviLimits = [.5 50]; end
                   
            tofSlider = uicontrol('Parent', [], 'Style', 'slider', ...
                'Value',  1, 'Min', log10(opt.pviLimits(1)), 'Max', log10(opt.pviLimits(2)), 'Visible', 'off');
            timeButton = uicontrol('Parent', [], 'Style', 'radiobutton', ...
                'Value',  0, 'Visible', 'off', 'Tag', 'years');
            timeText   = uicontrol('Parent', [], 'Style', 'text', ...
                'Value',  [], 'String', 'years', 'Visible', 'off');
            timeEdit   = uicontrol('Parent', [], 'Style', 'edit', ...
                'Value',  [], 'String', '', 'Visible', 'off', 'Tag', 'years');
            timePVButton = uicontrol('Parent', [], 'Style', 'radiobutton', ...
                'Value',  1, 'Visible', 'off', 'Tag', 'pvi');
            timePVText   = uicontrol('Parent', [], 'Style', 'text', ...
                'Value',  [], 'String', 'pvi', 'Visible', 'off');
            timePVEdit   = uicontrol('Parent', [], 'Style', 'edit', ...
                'Value',  [], 'String', num2str(opt.pviLimits(2)), 'Visible', 'off', 'Tag', 'pvi');
                                                              
            controls      = {{[], tofSlider, []}, {timeButton, timeText, timeEdit, timePVButton, timePVText, timePVEdit}};
            controlLayout = {[.1, .8, .1], [.06 .19 .25 .06 .19 .25]}; 
               
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                         'controlWidths', controlLayout, 'Title', opt.Title, ...
                         'Position', opt.Position, 'Visible', 'off', extraOpt{:});
            
            s.fixedHeight = true;
            % main callback
            [s.pviLimits, s.yearLimits] = deal(opt.pviLimits, opt.yearLimits);
            s.Callback = opt.Callback;
            s.tofSlider = tofSlider;
            s.yearEdit  = timeEdit;
            s.pviEdit   = timePVEdit;
            s.yearButton = timeButton;
            s.pviButton = timePVButton;
            % item callbacks (include main)
            tofSlider.Callback  = @s.sliderCallback;
            timeEdit.Callback   = @s.editCallback;
            timePVEdit.Callback = @s.editCallback;
            timeButton.Callback = @s.buttonCallback;
            timePVButton.Callback = @s.buttonCallback;
            %timeButton.Callback = @ 
            % set visible
            s.Visible = opt.Visible;
        end
        
        function val = get.Value(s)
            val = 10^s.tofSlider.Value;
        end
        function set.Value(s, val)
            if strcmp(s.timeUnit, 'pvi')
                s.pviEdit.String = num2str(val);
                s.editCallback(s.pviEdit, nan);
            else
                s.yearEdit.String = num2str(val);
                s.editCallback(s.yearEdit, nan);
            end
        end
            
        
        function val = get.timeUnit(s)
            if s.yearButton.Value == 1
                val = 'years';
            else
                val = 'pvi';
            end
        end
        
        function set.timeUnit(s,val)
            if ~strcmp(val, s.timeUnit)
                if strcmp(val, 'years')
                    s.yearButton.Value = 1;
                    s.pviButton.Value = 0;
                else
                    s.yearButton.Value = 0;
                    s.pviButton.Value = 1;
                end
                s.switchUnit;
            end
        end
        
        function sliderCallback(s, src, event)
            tof = 10^s.tofSlider.Value;
            if strcmp(s.timeUnit, 'years')
                s.yearEdit.String = sprintf('%0.4g', tof);
                s.pviEdit.String  = '';
            else
                s.pviEdit.String  = sprintf('%0.4g', tof);
                s.yearEdit.String = '';
            end
            s.Callback(src, event);
        end
        
        function editCallback(s, src, event)
            val = str2double(src.String);
            if isnan(val)
                s.sliderCallback(src, event);
            else
                if ~strcmp(s.timeUnit, src.Tag)
%                     if strcmp(src.Tag, 'years')
%                         s.pviButton.Value = ~src.Value;
%                     else
%                         s.yearButton.Value = ~src.Value;
%                     end
                    s.switchUnit(src.Tag);
                end
                s.tofSlider.Value = max(s.tofSlider.Min, min(s.tofSlider.Max, log10(val)));
                s.sliderCallback(src, event)
                s.Callback(src, event)
            end
        end
        
        function buttonCallback(s, src, event)
            %if strcmp(src.Tag, 'years')
            %    s.pviButton.Value = ~src.Value;
            %else
            %    s.yearButton.Value = ~src.Value;
            %end
            if (strcmp(src.Tag, 'years') && src.Value == 1) || ...
                    (strcmp(src.Tag, 'pvi') && src.Value == 0)
                unit = 'years';
            else
                unit = 'pvi';
            end
            s.switchUnit(unit);
            s.sliderCallback(src, event);
        end
            
            
        function switchUnit(s, unit)
            val = s.tofSlider.Value;
            if strcmp(unit, 'pvi')
                [lOld, lNew] = deal(s.yearLimits, s.pviLimits);
                s.pviButton.Value = 1;
                s.yearButton.Value = 0;
            else
                [lNew, lOld] = deal(s.yearLimits, s.pviLimits);
                 s.yearButton.Value = 1;
                 s.pviButton.Value = 0;
            end
            [lOld, lNew] = deal(log10(lOld), log10(lNew));
            w = (val-lOld(1))/diff(lOld);
            if ~(0<=w && w<=1)
                w =0;
            end
            newVal = (1-w)*lNew(1) + w*lNew(2); 
            set(s.tofSlider, {'Value', 'Min', 'Max'}, ...
                {newVal, lNew(1), lNew(2)});
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
