function moduleGUI()
%Interactive user interface for activation/deactivation of known mrst modules
%
% SYNOPSIS:
%  moduleGUI();
%  
% DESCRIPTION:
%   No input or output arguments. `moduleGUI` launches an interactive user
%   interface for activation and de-activation of models. In addition, some
%   rudimentary features for listing files by right-clicking a module is
%   included.
%
% SEE ALSO:
%   `mrstModule`, `mrstPath`

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

persistent f
% Make a figure with a button group for the modules, leaving some space for
% pushbuttons at the bottom.
if ~isempty(f) && ishandle(f)
    % Avoid creating more than one window, since we can't rely on close all
    % to get rid of them.
    close(f)
end

figNo = 100;
while ishandle(figNo)
    figNo = figNo + 1;
end
f = figure(figNo);
set(f,     'Name',          'MRST Module loader', ...
           'Toolbar',       'none', ...
           'NumberTitle',   'off', ...
           'MenuBar',       'none');
       
bg = uibuttongroup('Visible','on',...
                   'Position',[0 .1 1 .9]);

% Retrieve modules and sort them in alphabetic order
modules = sort(mrstPath());
Nc = numel(modules);

% Four rows in total for gui
nrow = 4;

wdth = 1/nrow;
hght = min(0.8*nrow/Nc, 0.8);

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

        m = modules{j};
        dims = [xoffset, yoffset, wdth, hght];
        hb(j) = uicontrol(bg, 'Style', 'checkbox', 'String', m,...
                           'Units', 'normalized', 'Position', dims, ...
                           'TooltipString', mrstPath('query', m), ...
                           'ButtonDownFcn', @buttonClicked, ...
                           'Callback', @buttonCallback);
    end

end
% Add button controls
addbutton = @(x, label, varargin) uicontrol(f, 'Units', 'normalized', ...
                                               'Position', [x, 0.01, 0.2, .08], ...
                                               'Style', 'pushbutton', ...
                                               'String', label, ...
                                               varargin{:});
% Button for unloading all modules
addbutton(0.05, 'Unload all', 'Callback', @clearModules)
% Show known modules in terminal
addbutton(0.2 + 3*0.025, 'List paths', 'Callback', @(src, event) mrstPath('list'))
% Update. Should go automatically for most cases, but since Matlab does not
% allow focus events officially, we include it.
addbutton(0.4 + 4*0.025, 'Update', 'Callback', @(src, event) updateButtons());
% Nice and large close button.
addbutton(0.6 + 5*0.025, 'Exit', 'Callback', @(src, event) close(f));

function updateButtons()
    % Update all buttons with currently active status
    active = getActive(modules);
    % Module number is low, use for loop for clarity and robustness
    for modNo = 1:Nc
        if ishandle(hb(modNo))
            set(hb(modNo), 'Value', active(modNo));
        end
    end
end

function buttonCallback(src, event) %#ok
    % Button has changed status
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

function buttonClicked(src, event)
    % Button was clicked (right mouse button or mousewheel).
    name = get(src, 'String');
    pth =  mrstPath(name);
    fprintf('* Directory listing for module %s \n', name)
    fprintf('* <a href="%s">%s</a>\n', pth, pth); 
    files = dir(pth);
    for ix = 1:numel(files)
        listStr(pth, files(ix))
    end
end

function clearModules(src, event) %#ok
    mrstModule clear
    updateButtons();
end
% Events to ensure that the buttons get updated when the user interacts
% with the window while having loaded modules through scripts/cmd
set(f, 'WindowKeyPressFc', @(src, event) updateButtons());
set(f, 'WindowButtonDownFcn', @(src, event) updateButtons());

updateButtons();
set(f, 'HandleVisibility',  'Callback')
end

function active = getActive(modules)
    current = mrstModule();
    active = cellfun(@(x) any(strcmpi(current, x)), modules);
end

function listStr(pth, x)
    name = x.name;
    if x.isdir
        if strncmp(name, '.', 1)
            return
        end
        fprintf('<a href="%s">%s</a>\n', fullfile(pth, name), name);
    else
        [tmp, tmp, ext] = fileparts(name); %#ok
        if strcmpi(ext, '.m')
            % Print edit link to .m-file
            fprintf('<a href="matlab:edit %s">%s</a>\n', fullfile(pth, name), name);
        else
            % Just print the name without a link, it isn't a m-file.
            fprintf('%s\n', name);
        end
    end
end
