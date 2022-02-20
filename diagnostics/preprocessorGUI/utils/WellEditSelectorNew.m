classdef WellEditSelectorNew < UIItem
    properties
        Callback
        W
        W_orig
        wellPopup
        typePopup
        controlPopup
        valueEdit
        valueText
        openCheck
        newButton
        resetButton
        saveButton
        newCount     = 0;
        prevSelection  = 1;
        changedCallback       
    end
    
    properties (Dependent)
        wellNo
        wellInfo
        sign
        type
        val
        status
        enableSaveReset % also used to detect changes
    end
    
    methods
        function w = WellEditSelectorNew(W, varargin)
            opt = struct('Parent',          [], ...
                         'Callback',        @(src, event)'', ...
                         'Position',        [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title',           'Edit/Add well', ...
                         'wellNames',       {{'w1', 'w2', 'w3'}}, ...
                         'controls',        {{'bhp', 'rate'}}, ...
                         'wellTypes',       {{'injector', 'producer'}}, ...
                         'changedCallback', @(src, event)'');
            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            getui = @(style, string, value, tag)uicontrol('Parent', [], 'Style', style, 'Value', value, ...
                                                  'Visible', 'off', 'String', string, 'Tag', tag);
            getuitext = @(string)getui('text', string, [], '');
            
            wellPopup    = getui('popup', opt.wellNames, 1, 'name');
            typeText     = getuitext('Type'); 
            typePopup    = getui('popup', opt.wellTypes, 1, 'sign');
            controlText  = getuitext('Control'); 
            controlPopup = getui('popup', opt.controls, 1, 'type');
            valueText    = getuitext('Value'); 
            valueEdit    = getui('edit', '100', [], 'val');
            openText     = getuitext('Open'); 
            openCheck    = getui('checkbox', '', 1, 'status');
            
%             newProdButton = getui('pushbutton', 'New producer', [], 'newprod');
%             newInjButton  = getui('pushbutton', 'New injector', [], 'newinj');
            
%             trajText        = getuitext('Trajectory:');                                   
            newButton   = getui('pushbutton', 'New', [], 'edit');
            resetButton = getui('pushbutton', 'Undo', [], 'reset');
            saveButton  = getui('pushbutton', 'Save', [], 'save');
            
            controls      = {{wellPopup}, ...
                             {typeText, typePopup}, ...
                             {controlText, controlPopup}, ...
                             {valueText, valueEdit, openText, openCheck}, ...
                             {newButton, resetButton, saveButton}};
            
            controlLayout = {1, ...
                             [.3, .7], ...
                             [.3, .7], ...
                             [.3, .3, .3, .1], ...
                             [.33, .34, .33]};
                         
           w = w@UIItem('Parent', opt.Parent, 'controls', controls, ...
                        'controlWidths', controlLayout, 'Title', opt.Title, ...
                        'Position', opt.Position, 'Visible', 'off', extraOpt{:});
           w.wellPopup       = wellPopup;
           w.typePopup       = typePopup;
           w.controlPopup    = controlPopup;
           w.valueEdit       = valueEdit;
           w.valueText       = valueText;
           w.openCheck  	 = openCheck;
           w.newButton       = newButton;
           w.resetButton     = resetButton;
           w.saveButton      = saveButton;
           
           w.Callback = @opt.Callback;
           w.wellPopup.Callback = @w.popupCallback;
           set([w.typePopup, w.controlPopup, w.valueEdit, w.openCheck], ...
               'Callback', @w.wellStructUpdateCallback)
           w.newButton.Callback   = @w.newCallback;
           w.resetButton.Callback = @w.resetCallback;
           w.saveButton.Callback  = @w.saveCallback;
           
           w.enableSaveReset  = 'off';
           if ~isempty(W)
               w.updateWells(W);
           end
           w.typePopup.Enable = 'off';
           % set visible
           w.Visible = opt.Visible;
        end
        
        function set.wellNo(w, val)
            w.wellPopup.Value = val;
        end
        function val = get.wellNo(w)
            val = w.wellPopup.Value;
        end
            
        function set.sign(w, val)
            if val == 1
                w.typePopup.Value = 1;
            elseif val == -1
                w.typePopup.Value = 2;
            else
                error('Well can''t have sign %d', val);
            end
            w.W(w.wellNo).sign = val;
        end
        function val = get.sign(w)
            if w.typePopup.Value == 1
                val = 1;
            else
                val = -1;
            end
        end
        
        function set.type(w, val)
            switch val
                case 'bhp'
                    w.controlPopup.Value = 1;
                    w.valueText.String = '[bar]';
                case 'rate'
                    w.controlPopup.Value = 2;
                    w.valueText.String = '[m3/day]';
                otherwise
                    error('Can''t set well-type %s', val);
            end
            w.W(w.wellNo).type = val;
            w.W(w.wellNo).val  = w.val;
        end
        function val = get.type(w)
            if w.controlPopup.Value == 1
                val = 'bhp';
            else
                val = 'rate';
            end
        end
        
        function set.val(w, val)
            switch w.type
                case 'bhp'
                    w.valueEdit.String = num2str(convertTo(val, barsa));
                case 'rate'
                    w.valueEdit.String = num2str(convertTo(val, 1/day));
            end
            w.W(w.wellNo).val = val;
        end
        function val = get.val(w)
            switch w.type
                case 'bhp'
                    val = convertFrom(str2double(w.valueEdit.String), barsa);
                case 'rate'
                    val = convertFrom(str2double(w.valueEdit.String), 1/day);
            end
        end
             
        function set.status(w, val)
            if islogical(val) && val
                w.openCheck.Value = 1;
            elseif islogical(val) && ~val
                w.openCheck.Value = 0;
            else
                error('Exspected logical, got %s', class(val));
            end
            w.W(w.wellNo).status = val;
        end
        function val = get.status(w)
            if w.openCheck.Value == 1
                val = true;
            else
                val = false;
            end
        end
        
        function set.enableSaveReset(w,val)
            if ~any(strcmp(val, {'on', 'off'}))
                error('Unkown enable option: %s', val);
            else
                w.saveButton.Enable  = val;
                w.resetButton.Enable = val;    
            end
            w.changedCallback([], []);
        end
        function val = get.enableSaveReset(w)
            if ~strcmp(w.saveButton.Enable, w.resetButton.Enable)
                % should be the same
                w.resetButton.Enable = w.saveButton.Enable;
            end
            val = w.saveButton.Enable;
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
            if strcmp(w.enableSaveReset, 'on')
                % reset
                w.wellNo = w.prevSelection;
                interuptMessage(w);
            else
                ix   = w.wellNo;
                well = w.W(ix);
                w.sign   = well.sign;
                w.type   = well.type;
                w.val    = well.val;
                w.status = well.status;
                w.Callback(src, event);
                w.prevSelection = ix;
            end
        end
        
        function wellStructUpdateCallback(w, src, event)
            switch src.Tag
                case 'sign'
                    w.sign = w.sign;
                case 'type'
                    w.type = w.type;
                case 'val'
                    w.val = w.val;
                case 'status'
                    w.status = w.status;
                otherwise
                    error('Unknown tag: %s', src.Tag)
            end
            w.enableSaveReset = 'on';
        end
        
        function newCallback(w, src, event)
            if strcmp(w.enableSaveReset, 'on')
                interuptMessage(w);
            else
                % copy current selection
                w.newCount = w.newCount + 1;
                w_new = w.W(w.wellNo);
                if w_new.sign > 0
                    nm = sprintf('INJ_NEW%d', w.newCount);
                else
                    nm = sprintf('PROD_NEW%d', w.newCount);
                end
                w_new.name = nm;
                w.W = vertcat(w.W, w_new);
                w.wellPopup.String = {w.W.name};
                w.enableSaveReset = 'on';
                % select new
                w.wellPopup.Value = numel(w.W);
                w.prevSelection   = numel(w.W);
                w.Callback(src, event);
            end
        end
        
        function resetCallback(w, src, event)
            dispif(mrstVerbose, 'Disregarding any changes since last save\n')
            w.enableSaveReset = 'off';
            if numel(w.W) > numel(w.W_orig)  % remove new well
                w.wellNo = 1;
                w.updateWells(w.W_orig);
                w.newCount = w.newCount - 1;
            else
                W_tmp = w.W;
                W_tmp(w.wellNo) = w.W_orig(w.wellNo);
                w.updateWells(W_tmp);
            end
            w.Callback(src, event);
        end
        
        function saveCallback(w, src, event)
            dispif(mrstVerbose, 'Saving current changes\n')
%             tmp = computeTraversedCells(G, traj, 'faces', opt.faces, ...
%                                     'exteriorFaceCorrection', true);
         
            w.W_orig = w.W;
            w.enableSaveReset = 'off';
            w.Callback(src, event);
        end
        
        function updateWells(w, W)
            if strcmp(w.enableSaveReset, 'on')
                interuptMessage(w);
            else
                w.W = W;
                % keep copy of un-edited
                w.W_orig = W;
                w.wellPopup.String = {W.name};
                w.popupCallback(w.wellPopup, []);
            end
        end
        
        function interuptMessage(w)
            fprintf('Well %s contains unsaved changes. Please save or undo before continuing\n', ...
                w.W(w.wellNo).name)
            col = w.saveButton.BackgroundColor;
            for k = 1:5
                set([w.saveButton, w.resetButton], 'BackgroundColor', mod(col - [.5 0 .5],1));
                pause(.1)
                set([w.saveButton, w.resetButton], 'BackgroundColor', col);
                pause(.1)
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
