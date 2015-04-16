function moduleGUI()
%Interactive user interface for activation/deactivation of known mrst modules
%
% SYNOPSIS:
%  moduleGUI();
%  
% PARAMETERS:
%   None
%
% RETURNS:
%   Nothing.
%
% SEE ALSO:
%   mrstModule, mrstPath

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


    f = figure('ToolBar', 'none');
    bg = uibuttongroup('Visible','on',...
                      'Position',[0 .1 1 .9]);

    % Put modules in alphabetic order
    modules = sort(mrstPath());

    Nc = numel(modules);

    % Four rows in total for gui
    nrow = 4;

    wdth = 1/nrow;
    hght = 0.1;

    perRow = ceil(Nc/nrow);
    dy = 1/perRow;
    
    % Place checkbox matrix for all modules
    hb = zeros(Nc, 1);
    for i = 1:nrow
        pp = 1;
        for j = ((i-1)*perRow + 1):min(Nc, i*perRow)
            xoffset = wdth*(i - 1);
            yoffset = 1 - pp*dy;

            pp = pp + 1;

            dims = [xoffset, yoffset, wdth, hght];
            hb(j) = uicontrol(bg, 'Style', 'checkbox', 'String', modules{j},...
                               'Units', 'normalized', 'Position', dims, ...
                               'Callback', @buttonPressed);
        end

    end

    addbutton = @(x, label, varargin) uicontrol(f, 'Units', 'normalized', ...
                                                   'Position', [x, 0.01, 0.2, .08], ...
                                                   'Style', 'pushbutton', ...
                                                   'String', label, ...
                                                   varargin{:});
    addbutton(0.05, 'Unload all', 'Callback', @clearModules)
    addbutton(0.2 + 3*0.025, 'List paths', 'Callback', @(src, event) mrstPath('list'))
    addbutton(0.4 + 4*0.025, 'Update', 'Callback', @(src, event) updateButtons());
    addbutton(0.6 + 5*0.025, 'Exit', 'Callback', @(src, event) close(f));

    function updateButtons()
        active = getActive(modules);
        for modNo = 1:Nc
            set(hb(modNo), 'Value', active(modNo));
        end
    end
    
    function buttonPressed(src, event)
        % Button pressed
        v = get(src, 'Value');
        name = get(src, 'String');
        if v
            % Value is positive, we add the module
            mrstModule('add', name);
            updateButtons();
        else
            % Value was negative, reset
            updateButtons()
            active = getActive(modules);
            active = active & ~strcmpi(modules, name);
            mrstModule('reset', modules{active});
            updateButtons();
        end
    end

    function clearModules(src, event)
        mrstModule clear
        updateButtons();
    end

    set(f, 'WindowKeyPressFc', @(src, event) updateButtons());
    set(f, 'WindowButtonDownFcn', @(src, event) updateButtons());
    
    updateButtons();
end

function active = getActive(modules)
    current = mrstModule();

    active = cellfun(@(x) any(strcmpi(current, x)), modules);
end

