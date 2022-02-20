classdef WellEditSelector < UIItem
    properties
        Callback
        W
        wellPopup
        %nameEdit
        typePopup
        controlPopup
        valueEdit
        newButton
        viewButton
        saveButton
        launchButton
    end
    properties (Dependent)
        wellNo
        wellInfo
    end
    
    methods
        function w = WellEditSelector(W, varargin)
            opt = struct('Parent',          [], ...
                         'Callback',        @(src, event)'', ...
                         'Position',        [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title',           'Edit/Add well', ...
                         'wellNames',       {{'w1', 'w2', 'w3'}}, ...
                         'controls',        {{'bhp', 'rate'}}, ...
                         'wellTypes',       {{'injector', 'producer'}});
                        % 'drawCallback',    @(src, event)disp('Boo'), ...
                        % 'newCallback',     @(src, event)disp('Boo'), ...
                        % 'viewCallback',  @(src, event)disp('Boo'), ...
                        % 'saveCallback)',   @(src, event)disp('Boo'));
                     [opt, extraOpt] = merge_options(opt, varargin{:});
            
            wellPopup = uicontrol('Parent', [], 'Style', 'popup', ...
                                  'Value',  1, 'String', opt.wellNames, 'Visible', 'off');
            
            %nameText  = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ...
            %                      'String', 'Name', 'Visible', 'off');
            %nameEdit  = uicontrol('Parent', [], 'Style', 'edit', ...
            %                      'String', 'sdfs', 'Visible', 'off');
            typeText  = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ...
                                  'String', 'Type', 'Visible', 'off');
            typePopup = uicontrol('Parent', [], 'Style', 'popup', ...
                                  'Value',  1, 'String', opt.wellTypes, 'Visible', 'off');
            controlText = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ...
                                    'String', 'Control', 'Visible', 'off');
            controlPopup = uicontrol('Parent', [], 'Style', 'popup', ...
                                     'Value',  1, 'String', opt.controls, 'Visible', 'off');
            valueText = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ...
                                  'String', 'Value ', 'Visible', 'off');
            valueEdit = uicontrol('Parent', [], 'Style', 'edit', ...
                                  'String', num2str(400), 'Visible', 'off');
                              
            newButton        = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'new', ...
                                         'String', 'New', 'Visible', 'off');
            viewButton     = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'view', ...
                                         'String', 'View', 'Visible', 'off');
            saveButton       = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'save', ...
                                         'String', 'Save', 'Visible', 'off');
            launchButton   = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'launch', ...
                                         'String', 'Launch diagnostics', 'Visible', 'off');
             
            %                 {nameText, nameEdit}, ...                        
            controls      = {{wellPopup}, ...
                             {typeText, typePopup}, ...
                             {controlText, controlPopup}, ...
                             {valueText, valueEdit}, ...
                             {newButton, viewButton, saveButton}, ...
                             {launchButton}};
            
            %                 [.4, .6], ...
            controlLayout = {1, ...
                             [.4, .6], ...
                             [.4, .6], ...
                             [.4, .6], ...
                             [.33, .34, .33], ...
                             1};
                         
            w = w@UIItem('Parent', opt.Parent, 'controls', controls, ...
                         'controlWidths', controlLayout, 'Title', opt.Title, ...
                         'Position', opt.Position, 'Visible', 'off', extraOpt{:});
           w.wellPopup = wellPopup;
           %w.nameEdit  = nameEdit;
           w.typePopup = typePopup;
           w.controlPopup = controlPopup;
           w.valueEdit    = valueEdit;
           w.newButton = newButton;
           w.viewButton = viewButton;
           w.saveButton = saveButton;
           w.launchButton = launchButton;
           
           w.Callback = @opt.Callback;
           w.wellPopup.Callback = @w.popupCallback;
           w.newButton.Callback  = @opt.Callback;
           w.viewButton.Callback = @opt.Callback;
           w.saveButton.Callback   = @opt.Callback;
           w.launchButton.Callback = @opt.Callback;
           
           if ~isempty(W)
               w.updateWells(W);
           end
           % set visible
           w.Visible = opt.Visible;
        end
        
        function set.wellNo(w, val)
            w.wellPopup.Value = val;
        end
        function val = get.wellNo(w)
            val = w.wellPopup.Value;
        end
        
        function s = get.wellInfo(w)
            nm = w.wellPopup.String{w.wellPopup.Value};
            tp = w.typePopup.String{w.typePopup.Value};
            if strcmp(tp, 'injector')
                sgn = 1;
            else
                sgn = -1;
            end
            cnt = w.controlPopup.String{w.controlPopup.Value};
            v = str2double(w.valueEdit.String);
            if strcmp(cnt, 'bhp')
                v = v*barsa;
            else
                v = v/day;
            end
            s = struct('name', nm, 'type', cnt, 'sign', sgn, 'val', v); 
        end
        
        function popupCallback(w, src, event)
            ix   = w.wellPopup.Value;
            well = w.W(ix);
            if well.sign > 1
                w.typePopup.Value = 1;
            else
                w.typePopup.Value = 2;
            end
            val = well.val;
            if strcmp(well.type, 'bhp')
                w.controlPopup.Value = 1;
                w.valueEdit.String = num2str(val/barsa);
            else
                w.controlPopup.Value = 2;
                w.valueEdit.String = num2str(abs(val)*day);
            end
            w.Callback(src, event);
        end
                
        function updateWells(w, W)
            w.W = W;
            w.wellPopup.String = {W.name};
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
