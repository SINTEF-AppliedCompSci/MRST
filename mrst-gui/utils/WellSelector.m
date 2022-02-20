classdef WellSelector < UIItem
    properties
        Callback
        injSelector
        prodSelector
        thresholdEdit   = [];
        thresholdSlider = [];
        autoCheckBox    = false;
        autoEdit        = [];
        communicationMatrix = [];
    end
    properties (Dependent)
        injectorIx
        producerIx
        threshold
        communicationLimit
        autoDetect
    end
    
    methods
        
        function s = WellSelector(varargin)
            
            opt = struct('Parent',          [], ...
                         'Callback',        'disp(''Hello'')', ...
                         'Position', [1 1 200 200], ...
                         'Visible',         'on', ...
                         'Title', 'Select interaction region', ...
                         'injectors', {{'I1', 'I2', 'I3'}}, ...
                         'producers', {{'P1', 'P2'}}, ...
                         'includeThreshold', true, ...
                         'includeTooltips', true, ...
                         'TooltipString', 'Right-click for options');
           [opt, extraOpt] = merge_options(opt, varargin{:});
           
            injHeader    = uicontrol('Parent', [], 'Style', 'text', ...
                                     'Value', [], 'String', 'Injectors', 'Visible', 'off');
            prodHeader   = uicontrol('Parent', [], 'Style', 'text', ...
                                     'Value', [], 'String', 'Producers', 'Visible', 'off');
            injSelector  = uicontrol('Parent', [], 'Style', 'listbox',  'Max', 2, 'Min', 0, ...
                                     'Value', [], 'String', opt.injectors, 'Visible', 'off');
            prodSelector = uicontrol('Parent', [], 'Style', 'listbox',  'Max', 2, 'Min', 0, ...
                                     'Value', [], 'String', opt.producers, 'Visible', 'off');
            autoCheckBox = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'Auto-detect well-pairs', ...
                                     'Visible', 'off');
            autoEdit     = uicontrol('Parent', [], 'Style', 'edit', ... 
                                     'String', num2str(1), 'Visible', 'off');
            if opt.includeTooltips
                injSelector.TooltipString  = opt.TooltipString;
                prodSelector.TooltipString = opt.TooltipString;
            end
            controls = { {injHeader,      prodHeader}, ...
                              {injSelector,    prodSelector}, ...
                              {autoCheckBox, autoEdit}};
            controlWidths = {[.5 .5], [nan nan], [nan nan]};
                          
            if opt.includeThreshold
                thresholdHeader = uicontrol('Parent', [], 'Style', 'text', ...
                                             'Value', [], 'String', 'Filter on lower threshold', ...
                                             'Visible', 'off');
                thresholdEdit = uicontrol('Parent', [], 'Style', 'edit', ...
                                          'String', num2str(0), 'Visible', 'off');
                thresholdSlider = uicontrol('Parent', [], 'Style', 'slider', ...
                                            'Min', 0, 'Max', 1, 'Value', 0, ...
                                            'SliderStep', [.01, .1], 'Visible', 'off');
                
                controls = [controls, ...
                           {{thresholdHeader}, ...
                           {thresholdEdit, thresholdSlider}}];
                controlWidths = [controlWidths, {nan, [nan nan]}];
            end
                
                                      
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, 'controlWidths', controlWidths, 'Title', opt.Title, ...
                         'Position', opt.Position,'Visible', 'off', extraOpt{:});
            
            
            % main callback
            s.Callback = opt.Callback;
            % item callbacks (include main)
            injSelector.Callback  = @s.injCallback;
            prodSelector.Callback = @s.prodCallback;
            autoCheckBox.Callback = @s.selectorCallback;
            autoEdit.Callback     = @s.autoCallback;
            
            s.injSelector     = injSelector;
            s.prodSelector    = prodSelector;
            s.autoCheckBox    = autoCheckBox;
            s.autoEdit        = autoEdit;
            s.communicationLimit = str2double(autoEdit.String);

            
            if opt.includeThreshold
                thresholdEdit.Callback = @s.editCallback;
                thresholdSlider.Callback = @s.sliderCallback;
                
                s.thresholdEdit   = thresholdEdit;
                s.thresholdSlider = thresholdSlider;
            end
            % item context menues
            
            
            % set visible
            s.Visible = 'on';
            s.Enable  = 'off';
            
            % add item context menues after figure has been created
            injSelector.UIContextMenu  = s.listboxContextMenu({'Select all injectors', 'Clear selected injectors'}, ...
                                                             {@s.selectAllInjectors, @s.clearAllInjectors} );
            prodSelector.UIContextMenu = s.listboxContextMenu({'Select all producers', 'Clear selected producers'}, ...
                                                             {@s.selectAllProducers, @s.clearAllProducers} );
 
        end
        
        function com = getCommunication(s)
           com = s.communicationMatrix;
           tot = sum(com(:));
           n   = max(size(com));
           com = com > (s.communicationLimit/100)*tot/n;
        end
        
        function cstr = getCommunicationStrength(s)
           cstr = s.communicationMatrix;
           cstr(~s.getCommunication()) = 0;
        end
        
        function prodIx = getConnectedProd(s,injIx)
           com = s.getCommunication();
           ix = logical(sum(com(injIx,:),1));
           prodIx = find(ix);
        end
        
        function injIx = getConnectedInj(s, prodIx)
           com = s.getCommunication();
           ix = logical(sum(com(:,prodIx), 2));
           injIx = find(ix).';
        end
        
        function set.injectorIx(s, val)
            s.injSelector.Value = val;
        end
        function val = get.injectorIx(s)
            val = s.injSelector.Value;
        end
        
        function set.producerIx(s, val)
            s.prodSelector.Value = val;
        end
        function val = get.producerIx(s)
            val = s.prodSelector.Value;
        end
        
        function set.threshold(s, val)
            s.thresholdSlider.Value = val;
            s.thresholdEdit.String  = num2str(val,'%7.2g');
        end
        function val = get.threshold(s)
            val = s.thresholdSlider.Value;
        end
            
        function set.autoDetect(s, val)
            s.autoCheckBox.Value = val;
        end
        function val = get.autoDetect(s)
            val = s.autoCheckBox.Value;
        end
        
        function set.communicationLimit(s, val)
           val = capValue(val, [0 100]);
           if isfinite(val)
               s.autoEdit.String = num2str(val);
           end
        end
        function val = get.communicationLimit(s)
           val = str2double(s.autoEdit.String);
        end
           
        function selectAllInjectors(s, src, event)
            nv = numel(s.injSelector.String);
            if numel(s.injectorIx) ~= nv
                s.injectorIx = (1:nv);
                s.injCallback(src, event)
            end
        end
        function selectAllProducers(s, src, event)
            nv = numel(s.prodSelector.String);
            if numel(s.producerIx) ~= nv
                s.producerIx = (1:nv);
                s.prodCallback(src, event)
            end
        end
        
        function clearAllInjectors(s, src, event)
            if numel(s.injectorIx) ~= 0
                s.injectorIx = [];
                if s.autoDetect && ~isempty(s.communicationMatrix)
                    s.producerIx = [];
                end
                s.selectorCallback(src, event);
            end
        end
        
        function clearAllProducers(s, src, event)
            if numel(s.producerIx) ~= 0
                s.producerIx = [];
                if s.autoDetect && ~isempty(s.communicationMatrix)
                    s.injectorIx = [];
                end
                s.selectorCallback(src, event);
            end
        end
        
        function editCallback(s, src, event)
            if all(ismember(s.thresholdEdit.String, '0123456789-+eE.'))
                s.thresholdSlider.Value = max(0, min(1, str2double(s.thresholdEdit.String)));
                s.selectorCallback(src, event);
            end
            s.threshold = s.thresholdSlider.Value;
        end
        
        function autoCallback(s, src, event)
           if all(ismember(s.autoEdit.String, '0123456789-+eE.'))
              s.autoEdit.Value = max(0, min(100,str2double(s.autoEdit.String)));
              s.selectorCallback(src,event);
           end
        end
        
        function sliderCallback(s, src, event)
            s.threshold = s.thresholdSlider.Value;
            s.selectorCallback(src, event);
        end
        
        function injCallback(s, src, event)
            if s.autoDetect && ~isempty(s.communicationMatrix)
                s.producerIx = s.getConnectedProd(s.injectorIx);
            end
            s.selectorCallback(src, event);
        end
        
        function prodCallback(s, src, event)
            if s.autoDetect && ~isempty(s.communicationMatrix)
                s.injectorIx = s.getConnectedInj(s.producerIx);
            end
            s.selectorCallback(src, event);
        end
        
        function selectorCallback(s, src, event)
            %main callback
            % XXXX inclide if statement to check for actual change
            if ischar(s.Callback)
                eval(s.Callback);
            else
                s.Callback(src, event);
            end
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
