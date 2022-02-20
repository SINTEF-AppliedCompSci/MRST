classdef TwoPhaseFluidSelector < UIItem
    properties
        Callback
        origButton
        customButton
        viscWaterEdit
        viscOilEdit
        satWaterEdit
        satOilEdit
        exponentWatEdit
        exponentOilEdit
        updateButton
        isTouched = true;
        %Axes
    end
    properties (Dependent)
        useOriginal
        muw
        muo
        swc
        soc
        n
    end
    
    methods
        
        function s = TwoPhaseFluidSelector(varargin)
            opt = struct('Parent',          [], ...
                         'Callback',        @(src, event)disp('Boo'), ...
                         'Position',        [10 10 300 200], ...
                         'Visible',         'on', ...
                         'Title',           'Two-phase fluid options', ...
                         'muw',             .3*centi*poise, ...
                         'muo',              1*centi*poise, ...
                         'swc',              0, ...
                         'soc',              0, ...
                         'n',                2);
            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            textOpts = {'Parent', [], 'Style', 'text',    'Value',  [], 'Visible', 'off'};
            editOpts = {'Parent', [], 'Style', 'edit',    'Value',  [], 'Visible', 'off'};
            
            origButton = uicontrol('Parent', [], 'Style', 'radiobutton', ...
                                   'Value',  1, 'Visible', 'off', 'Tag', 'orig');
            origText   = uicontrol('Parent', [], 'Style', 'text', ...
                                   'Value',  [], 'String', 'Use original', 'Visible', 'off');
            customButton = uicontrol('Parent', [], 'Style', 'radiobutton', ...
                                   'Value',  0, 'Visible', 'off', 'Tag', 'custom');
            customText   = uicontrol('Parent', [], 'Style', 'text', ...
                                   'Value',  [], 'String', 'Use custom', 'Visible', 'off');
            
            viscText   = uicontrol(textOpts{:}, 'String', 'Viscosities [cp]:');
            viscwText   = uicontrol(textOpts{:}, 'String', 'W');
            viscWaterEdit = uicontrol(editOpts{:}, 'String', num2str(opt.muw/(centi*poise)), 'Enable', 'off');
            viscoText   = uicontrol(textOpts{:}, 'String', 'O');
            viscOilEdit = uicontrol(editOpts{:}, 'String', num2str(opt.muo/(centi*poise)), 'Enable', 'off');
            
            satText   = uicontrol(textOpts{:}, 'String', 'Critical sat.:');
            satwText   = uicontrol(textOpts{:}, 'String', 'W');
            satWaterEdit = uicontrol(editOpts{:}, 'String', num2str(opt.swc), 'Enable', 'off');
            satoText   = uicontrol(textOpts{:}, 'String', 'O');
            satOilEdit = uicontrol(editOpts{:}, 'String', num2str(opt.soc), 'Enable', 'off');
            
            expText    = uicontrol(textOpts{:}, 'String', 'Rel-perm exponent:');
            expWText   = uicontrol(textOpts{:}, 'String', 'W');
            exponentWatEdit = uicontrol(editOpts{:}, 'String', num2str(opt.n), 'Enable', 'off');
            expOText        = uicontrol(textOpts{:}, 'String', 'O');
            exponentOilEdit = uicontrol(editOpts{:}, 'String', num2str(opt.n), 'Enable', 'off');
            
            
            updateButton = uicontrol('Parent', [], 'Style', 'pushbutton', 'String', 'Update', ...
                                     'Visible', 'off', 'Enable', 'off');         

                               
            controls      = {{origButton, origText, customButton, customText}, ...
                             {viscText, viscwText, viscWaterEdit, viscoText, viscOilEdit}, ...
                             {satText, satwText, satWaterEdit, satoText, satOilEdit}, ...
                             {expText, expWText, exponentWatEdit, expOText,  exponentOilEdit}, ...
                             {[], updateButton}};
            controlLayout = {[.06, .44, .06, .44], ...
                             [.5, .1, .15, .1, .15], ...
                             [.5, .1, .15, .1, .15], ...
                             [.5, .1, .15, .1, .15], ...
                             [.7, .3]}; 
      
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                                'controlWidths', controlLayout, 'Title', opt.Title, ...
                                'Position', opt.Position,'Visible', 'off', extraOpt{:});
            
            s.origButton = origButton;
            s.customButton = customButton;
            s.viscWaterEdit = viscWaterEdit;
            s.viscOilEdit   = viscOilEdit;
            s.satWaterEdit  = satWaterEdit;
            s.satOilEdit    = satOilEdit;
            s.exponentWatEdit  = exponentWatEdit;
            s.exponentOilEdit  = exponentOilEdit;
            s.updateButton  = updateButton;        
            s.fixedHeight = true;
            
            % main callback
            s.Callback = opt.Callback;
            % item callbacks (include main)
            s.origButton.Callback    = @s.buttonCallback;
            s.customButton.Callback  = @s.buttonCallback;
            s.updateButton.Callback  = @s.updateCallback;
            s.viscWaterEdit.Callback = @(src, event)s.editCallback(src, event, 'muw');
            s.viscOilEdit.Callback   = @(src, event)s.editCallback(src, event, 'muo');
            s.satWaterEdit.Callback  = @(src, event)s.editCallback(src, event, 'swc');
            s.satOilEdit.Callback	 = @(src, event)s.editCallback(src, event, 'soc');
            s.exponentWatEdit.Callback  = @(src, event)s.editCallback(src, event, 'nw');
            s.exponentOilEdit.Callback  = @(src, event)s.editCallback(src, event, 'no');
            % set visible
            s.Visible = opt.Visible;
        end
        
        function val = get.useOriginal(s)
            val = (s.origButton.Value == 1);
        end
        
        function set.useOriginal(s, val)
            s.origButton.Value = val;
            s.buttonCallback(s.origButton, [])
        end
        
        function set.muw(s, val)
            val = capValue(val, [.01, 500]);
            if isfinite(val)
                s.viscWaterEdit.String = num2str(val);
            end
        end
        function val = get.muw(s)
            val = str2double(s.viscWaterEdit.String);
        end
        
        function set.muo(s, val)
            val = capValue(val, [.01, 500]);
            if isfinite(val)
                s.viscOilEdit.String = num2str(val);
            end
        end
        function val = get.muo(s)
            val = str2double(s.viscOilEdit.String);
        end
        
        function set.swc(s, val)
            val = capValue(val, [0, .9]);
            if isfinite(val)
                val = capValue(min(val, 1-s.soc-.1), [0 1]);
                s.satWaterEdit.String = num2str(val);
            end
        end
        function val = get.swc(s)
            val = str2double(s.satWaterEdit.String);
        end
        
        function set.soc(s, val)
            val = capValue(val, [0, .9]);
            if isfinite(val)
                val = capValue(min(val, 1-s.swc-.1), [0 1]);
                s.satOilEdit.String = num2str(val);
            end
        end
        function val = get.soc(s)
            val = str2double(s.satOilEdit.String);
        end
        
        function set.n(s, val)
            val = capValue(val, [1, 5]);
            if ~isempty(val)
                s.exponentWatEdit.String = num2str(val(1));
                s.exponentOilEdit.String = num2str(val(2));
            end
        end
        function val = get.n(s)
            val = [str2double(s.exponentWatEdit.String), ...
                   str2double(s.exponentOilEdit.String)];
        end
        
        function editCallback(s, src, event, fld)
            % trigger set-function to deal with invalid input
            v = str2double(src.String);
            if strcmp(fld, 'nw')
                s.n(1) = v;
            elseif strcmp(fld, 'no')
                s.n(2) = v;
            else
                s.(fld) = v;
            end
            s.isTouched = true;
        end
        
        function updateCallback(s, src, event)
            s.Callback(src, event)
            s.isTouched = false;
        end
        
        function buttonCallback(s, src, event)
            if (strcmp(src.Tag, 'orig') && src.Value == 1) || ...
               (strcmp(src.Tag, 'custom') && src.Value == 0)
                st = 'off';
                s.origButton.Value   = 1;
                s.customButton.Value = 0;
            else
                st = 'on';
                s.origButton.Value   = 0;
                s.customButton.Value = 1;
            end
            
            uic = {s.viscWaterEdit, s.viscOilEdit, s.satWaterEdit, ...
                   s.satOilEdit, s.exponentWatEdit, s.exponentOilEdit};
            for k = 1:numel(uic)
                uic{k}.Enable = st;
            end
            s.isTouched = true;
        end
        
    end
end

function v = capValue(v, lims)
if ischar(v)
    v = str2double(v);
end
assert(isnumeric(v), 'Non-numeric value...')
v = min(lims(2), max(lims(1), v));
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
