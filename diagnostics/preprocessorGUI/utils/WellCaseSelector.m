classdef WellCaseSelector < UIItem
    properties
        Callback
        listbox
        ixMax=0
    end
    properties (Dependent)
        ix
        selection
        names
    end
    
    methods
        
        function s = WellCaseSelector(varargin)
            sampleProps = {{{'base'}}};
            opt = struct('Parent',          [], ...
                         'Callback',        @(src, event)disp('Boo'), ...
                         'Position',        [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title',           'Well cases', ...
                         'names',           sampleProps);

            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            listbox = uicontrol('Parent', [], 'Style', 'listbox', 'Min', 1, 'Max', 10, ...
                              'Value',  1, 'String', opt.names{:}, 'Visible', 'off', 'Tag', 'select');
            newButton      = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'new', ...
                                         'String', 'New', 'Visible', 'off');
            deleteButton   = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'delete', ...
                                         'String', 'Delete', 'Visible', 'off');
            launchButton   = uicontrol('Parent', [], 'Style', 'pushbutton', 'Tag', 'launch', ...
                                         'String', 'Launch diagnostics', 'Visible', 'off');
            controls      = {{listbox}, ...
                             {newButton, deleteButton, launchButton}};
            controlLayout = {1, [.25 .25 .5]}; 
               
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                         'controlWidths', controlLayout, 'Title', opt.Title, ...
                         'Position', opt.Position, 'Visible', 'off', extraOpt{:});
            
            s.names = opt.names;
            s.listbox      = listbox;
            s.fixedHeight = true;
            % main callback
            s.Callback = opt.Callback;
            % item callbacks (include main)
            listbox.Callback        = @s.Callback;
            newButton.Callback    = @s.newCallback;
            deleteButton.Callback = @s.deleteCallback;
            launchButton.Callback  = @s.Callback;
            % set visible
            s.Visible = opt.Visible;
        end
        
        function set.ix(s, val)
            s.listbox.Value = val;
        end
        function val = get.ix(s)
            val = s.listbox.Value;
        end
        
        function val = get.selection(s)
            val = s.listbox.String{s.ix};
        end
        
        function set.names(s, val)
            s.listbox.String = val;
        end
        function val = get.names(s)
            val = s.listbox.String;
        end
        
        function newCallback(s, src, event)
            nn = numel(s.names);
            if nargin < 4
                s.ixMax = s.ixMax+1;
                name    = sprintf('case_%d', s.ixMax);
            end
            s.names = [s.names; {name}];
            s.Callback(src, event);
            s.ix = nn+1;
        end
        
        function deleteCallback(s, src, event)
            if any(s.ix == 1)
                warning('Can''t delete case %s', s.names{1});
            else
                tmp = s.ix;
                s.Callback(src, event);
                if tmp==numel(s.names)
                    s.ixMax=s.ixMax-1;
                end
                s.names(tmp) = [];
                s.ix = min(tmp)-1;
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
