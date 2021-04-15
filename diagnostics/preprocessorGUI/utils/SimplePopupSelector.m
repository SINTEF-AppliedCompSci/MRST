classdef SimplePopupSelector < UIItem
    properties
        Callback
        popup
        props
    end
    properties (Dependent)
        ix
        selection
    end
    
    methods
        
        function s = SimplePopupSelector(varargin)
            
            sampleProps = {{{'none',  'v2', 'v2'}}};
            opt = struct('Parent',          [], ...
                         'Callback',        @(src, event)disp('Boo'), ...
                         'Position',        [10 10 200 200], ...
                         'Visible',         'on', ...
                         'Title',           'Line plot', ...
                         'props',           sampleProps);

            [opt, extraOpt] = merge_options(opt, varargin{:});
            
            popup = uicontrol('Parent', [], 'Style', 'popup', ...
                              'Value',  1, 'String', opt.props{:}, 'Visible', 'off');
                                                              
            controls      = {{popup, []}};
            controlLayout = {[nan, .2]}; 
               
            s = s@UIItem('Parent', opt.Parent, 'controls', controls, ...
                         'controlWidths', controlLayout, 'Title', opt.Title, ...
                         'Position', opt.Position, 'Visible', 'off', extraOpt{:});
            
            s.props = opt.props;                            
            s.popup      = popup;
            s.fixedHeight = true;
            % main callback
            s.Callback = opt.Callback;
            % item callbacks (include main)
            popup.Callback      = @s.Callback;

            % set visible
            s.Visible = opt.Visible;
        end
        
        function set.ix(s, val)
            s.popup.Value = val;
        end
        function val = get.ix(s)
            val = s.popup.Value;
        end
        
        function val = get.selection(s)
            val = s.popup.String{s.ix};
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
