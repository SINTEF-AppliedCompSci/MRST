function varargout = mrstExploreModules()
% Interactively explore MRST modules and corresponding examples
%
% SYNOPSIS:
%   mrstExploreModules();
%
% DESCRIPTION:
%   Launches an interactive graphical user interface for examining MRST
%   modules, their examples and relevant papers.
%
% REQUIRED PARAMETERS:
%   None.
%
% RETURNS:
%   h  - Figure handle for the window (if requested).
%
% SEE ALSO:
%   mrstReferencesGUI, moduleGUI, mrstModule, mrstExamples

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
    df = get(0, 'DefaultFigurePosition');
    
    modules = sort(mrstPath());
    modules = ['MRST'; modules];
    lwidth = 0.20;
    midwidth = 1 - 2*lwidth;
    midpos = [lwidth, 0, midwidth, .5];
    [examples, examplenames] = deal([]);
    f = figure;
    set(f,     'Name',          'MRST module explorer', ...
               'Toolbar',       'none', ...
               'NumberTitle',   'off', ...
               'Position',      df.*[1, 1, 2, 1], ...
               'Units', 'normalized', ...
               'MenuBar',       'none');
    % Left column: List of modules
    modlist = uicontrol(f, 'Units', 'normalized', ...
                 'Position', [0, 0, lwidth, 1],...
                 'Style', 'listbox',...
                 'String', modules, ...
                 'Callback', @onClickModule);
    % Middle column: Title and description
    titl = uicontrol(f, 'Units', 'normalized', ...
                 'Position', [midpos(1), .9, midpos(3), .1],...
                 'Style', 'text',...
                 'min', 0, 'max', inf, ...
                 'HorizontalAlignment', 'center', 'FontSize', 20);

    descr = uicontrol(f, 'Units', 'normalized', ...
                 'Position', [midpos(1), .6, midpos(3), .325],...
                 'Style', 'edit',...
                 'min', 0, 'max', inf, ...
                 'HorizontalAlignment', 'left');
    % Reference GUI
    uicontrol(f, 'Units', 'normalized', ...
                 'Position', [midpos(1), .5, midpos(3), .075],...
                 'Style', 'text',...
                 'String', 'Relevant papers', ...
                 'min', 0, 'max', inf, ...
                 'HorizontalAlignment', 'center', 'FontSize', 16);
    paperh = mrstReferencesGUI([], f, midpos);
    
    % Right column: List of examples and corresponding buttons
    exlist = uicontrol(f, 'Units', 'normalized', ...
                 'Position', [midwidth + lwidth, .1, lwidth, .9],...
                 'Style', 'listbox');

    uicontrol(f, 'Units', 'normalized', ...
                 'Position', [midwidth + lwidth, 0, lwidth/2, .1],...
                 'Style', 'pushbutton',...
                 'String', 'Edit example', ...
                 'Callback', @onClickExample);
             
    uicontrol(f, 'Units', 'normalized', ...
                 'Position', [midwidth + 3/2*lwidth, 0, lwidth/2, .1],...
                 'Style', 'pushbutton',...
                 'Tag',    'publishbutton', ...
                 'String', 'Publish example', ...
                 'Callback', @onClickExample);

    function onClickModule(src, event)
        ix = get(src, 'Value');
        mods = get(src, 'String');
        module = mods{ix};
        if strcmpi(module, 'mrst')
            % Core mode
            m = [];
            module = 'core';
            name = 'MRST Core';
        else
            set(descr, 'String', getDescription(module))
            m = {module};
            name = ['Module "', module, '"'];
        end
        delete(paperh);
        set(titl, 'String', name);
        
        [examples, examplenames] = getExamples(module);
        set(exlist, 'String', examplenames);
        set(exlist, 'Value', 1);
        paperh = mrstReferencesGUI(m, f, midpos);
    end

    function onClickExample(src, event)
        ix = get(exlist, 'Value');
        if isempty(examples)
            return
        end
        fn = examples{ix};
        if strcmpi(get(src, 'Tag'), 'publishbutton')
            modNo = get(modlist, 'Value');
            mrstModule('add', modules{modNo});
            webfile = publish(fn);
            web(webfile);
        else
            edit(fn);
        end
    end

    onClickModule(modlist, []);
    
    if nargout > 0
        varargout{1} = f;
    end
end

function d = getDescription(module)
    if strcmpi(module, 'mrst')
        pth = ROOTDIR();
    else
        pth = mrstPath(module);
    end
    files = dir(pth);
    d = 'No description available.';
    for i = 1:numel(files)
        if strcmpi(files(i).name, 'readme') || strcmpi(files(i).name, 'readme.txt')
            d = fileread(fullfile(pth, files(i).name));
            break;
        end
    end
end

function [examples, names] = getExamples(module)
    ex = mrstExamples(module);
    
    examples = ex{1};
    nex = numel(examples);
    names = cell(nex, 1);
    for i = 1:nex
        spl = strsplit(examples{i}, filesep);
        sub = find(strcmpi(spl, 'examples')) + 1;
        names{i} = fullfile(spl{sub:end});
    end
end

