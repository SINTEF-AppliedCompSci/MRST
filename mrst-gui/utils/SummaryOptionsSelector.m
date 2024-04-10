classdef SummaryOptionsSelector < UIItem
    properties
        caseEdits
        caseSwitches
        casePropType
        namePropType
        propPropType
        subItemSwitch
        scaleSwitch
    end
    properties (Dependent)
        caseNames
        caseSelection
        doScale
        doSelectSubItems
        lineProperties
        Callback
    end
    
    methods 
        function s = SummaryOptionsSelector(varargin)
            opt = struct('Parent',          [], ...
                         'Callback',        @(varargin)disp('Hello!'), ...
                         'Position', [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title', 'Options', ...
                         'caseNames', {{'case1', 'case2'}});
            [opt, extraOpt] = merge_options(opt, varargin{:});
            n = numel(opt.caseNames);
            [caseEdits, caseSwitches] = deal(cell(1,n));
            for k = 1:n
                caseEdits{k}    = uicontrol('Parent', [], 'Style', 'edit', 'Visible', 'off', ...
                                      'String', opt.caseNames{k}, 'Min', 0, 'Max', 1);
                caseSwitches{k} = uicontrol('Parent', [], 'Style', 'checkbox', 'String', '', ...
                                         'Visible', 'off', 'Value', 1);
            end
            subItemSwitch = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'Sub-items', ...
                                         'Visible', 'off', 'Value', 1);
            scaleSwitch   = uicontrol('Parent', [], 'Style', 'checkbox', 'String', 'Scale', ...
                                         'Visible', 'off');
            caseText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                  'String', 'Cases:', 'Visible', 'off');
            nameText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                  'String', 'Names:', 'Visible', 'off');
            propText   = uicontrol('Parent', [], 'Style', 'text', 'Value', [], ... 
                                  'String', 'Props:', 'Visible', 'off');                              
            casePropType = uicontrol('Parent', [], 'Style', 'popup', 'String', {'Marker', 'Color', 'LineStyle'}, ...
                                     'Value', 1, 'Visible', 'off');
            namePropType = uicontrol('Parent', [], 'Style', 'popup', 'String', {'Marker', 'Color', 'LineStyle'}, ...
                                    'Value', 3, 'Visible', 'off');
            propPropType = uicontrol('Parent', [], 'Style', 'popup', 'String', {'Marker', 'Color', 'LineStyle'}, ...
                                     'Value', 2, 'Visible', 'off');  
            c1 = applyFunction(@(x,y){x,y}, caseEdits, caseSwitches);
            controls      = [c1, {{subItemSwitch, scaleSwitch}, {caseText,nameText, propText}, ...
                             {casePropType, namePropType, propPropType}}];
            cl1 = repmat({[.9, nan]}, [1, numel(caseEdits)]);
            controlLayout = [cl1, {[.5 nan], [.33, .33, nan], [.33, .33, nan]}];
                
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                         'controlWidths', controlLayout, 'Title', opt.Title, ...
                         'Position', opt.Position,'Visible', 'off', extraOpt{:});
            
            s.caseEdits = caseEdits;
            s.caseSwitches = caseSwitches;
            s.casePropType = casePropType;
            s.namePropType = namePropType;
            s.propPropType = propPropType;
            s.subItemSwitch = subItemSwitch;
            s.scaleSwitch   = scaleSwitch;
            
            s.Callback = opt.Callback;
            setCallbackForAll(s.controls, opt.Callback);
            
            s.Visible = opt.Visible;
        end
        
        function set.Callback(s, fn)
            setCallbackForAll(s.controls, fn);
        end
        
        function fn = get.Callback(s)
            fn = s.caseEdits{1}.Callback;
        end
        
        function set.caseNames(s, nameList)
            assert(numel(nameList) == numel(s.caseEdits));
            for k = 1:numel(nameList)
                s.caseEdits{k}.String = nameList{k};
            end
        end
        
        function val = get.caseNames(s)
            val = applyFunction(@(x)x.String, s.caseEdits);
        end
        
        function set.caseSelection(s, val)
            for k = 1:numel(val)
                s.caseSwitches{k}.Value = double(val(k));
            end
        end
        
        function val = get.caseSelection(s)
            val = logical(cellfun(@(x)x.Value, s.caseSwitches));
        end
       
        function set.lineProperties(s, typeList)
            tp    = {'casePropType', 'namePropType', 'propPropType'};
            for k = 1:numel(typeList)
                ii = strcmp(typeList{k}, s.(tp{k}).String);
                if nnz(ii) == 1
                    s.(tp{k}).Value = find(ii);
                end
            end
        end
        
        function val = get.lineProperties(s)
            tp  = {'casePropType', 'namePropType', 'propPropType'};
            ii  = applyFunction(@(tp)s.(tp).Value, tp);
            val = applyFunction(@(k, tp)s.(tp).String{k}, ii, tp);
        end
        
        function set.doScale(s, val)
            s.scaleSwitch.Value = double(val);
        end
        
        function val = get.doScale(s)
            val = logical(s.scaleSwitch.Value);
        end
        
        function set.doSelectSubItems(s, val)
            s.subItemSwitch.Value = double(val);
        end
        
        function val = get.doSelectSubItems(s)
            val = logical(s.subItemSwitch.Value);
        end
    end
end

function setCallbackForAll(c, fn)
if iscell(c)
    for k = 1:numel(c)
        setCallbackForAll(c{k}, fn);
    end
elseif isa(c, 'matlab.ui.control.UIControl')
    c.Callback = fn;
end
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
