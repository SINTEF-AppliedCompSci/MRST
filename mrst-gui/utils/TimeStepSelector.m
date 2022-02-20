classdef TimeStepSelector < UIItem
    properties
        Callback
        selector
        linkCallback = '';
    end
    properties (Dependent, SetObservable)
        Value = [];
        ix
    end
    
    methods
        
        function s = TimeStepSelector(varargin)
            
            opt = struct('Parent',          [], ...
                         'Callback',        'disp(''Hello'')', ...
                         'Position', [1 1 200 200], ...
                         'Visible',         'on', ...
                         'Title', 'Select time-steps', ...
                         'tSteps', {{'t1', 't2', 't3'}}, ...
                         'includeTooltips', true, ...
                         'TooltipString', 'Right-click for options');
            [opt, extraOpt] = merge_options(opt, varargin{:});
                                   
            selector    = uicontrol('Parent', [], 'Style', 'listbox',  'Max', 2, 'Min', 0, ...
                                    'Value', [], 'String', opt.tSteps, 'Visible', 'off');
            if opt.includeTooltips
                selector.TooltipString = opt.TooltipString;
            end
            
            s = s@UIItem('Parent', opt.Parent, 'controls', {{selector}}, 'Title', opt.Title, ...
                               'Position', opt.Position,'Visible', 'off', extraOpt{:});
            
            s.selector = selector;
            % main callback
            s.Callback = opt.Callback;
            selector.Callback  = @s.selectorCallback;
            
            % set visible
            s.Visible = opt.Visible;
            
            % add item context menues after figure has been created
            s.selector.UIContextMenu  = s.listboxContextMenu({'Select all steps', 'Clear selected steps'}, ...
                                                               {@s.selectAll, @s.clearAll} ); 
        end
        
        function set.Value(s, val)
            nv = numel(s.selector.String);
            assert(all(val>0) && all(val<=nv));
            s.selector.Value = val;
        end
        function val = get.Value(s)
            val = s.selector.Value;
        end
        
        function set.ix(s, val)
            s.Value = val;
        end
        function val = get.ix(s)
            val = s.Value;
        end
        
        function selectAll(s, src, event)
            nv = numel(s.selector.String);
            if nv ~= numel(s.Value)
                s.Value = (1:nv);
                s.selectorCallback(src, event)
            end
        end
        
        function clearAll(s, src, event)
            if ~isempty(s.Value)
                s.Value = [];
                s.selectorCallback(src, event)
            end
        end
        
        function selectorCallback(s, src, event)
            %main callback
            if ~isempty(s.linkCallback)
                s.linkCallback(src, event)
            end
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
