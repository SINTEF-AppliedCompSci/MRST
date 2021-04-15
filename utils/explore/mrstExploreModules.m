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
% RETURNS:
%   h  - Figure handle for the window (if requested).
%
% SEE ALSO:
%   `mrstReferencesGUI`, `moduleGUI`, `mrstModule`, `mrstExamples`

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
    if ~isempty(f) && ishandle(f)
        % Avoid creating more than one window, since we can't rely on close all
        % to get rid of them.
        close(f)
    end

    figNo = 100;
    while ishandle(figNo)
        figNo = figNo + 1;
    end
    modules = sort(mrstPath());
    modules = ['MRST'; modules];
    lwidth = 0.20;
    midwidth = 1 - 2*lwidth;
    topfrac = 0.66;
    topPos = [lwidth, 1 - topfrac, midwidth, topfrac];
    endPos = [lwidth, 0, midwidth, 1 - topfrac];
    
    [examples, examplenames] = deal([]);
    figpos = getPanelPosition();
    f = figure(figNo);
    set(f,     'Name',          'MRST module explorer', ...
               'Toolbar',       'none', ...
               'NumberTitle',   'off', ...
               'Position',      figpos, ...
               'Units',         'normalized', ...
               'MenuBar',       'none');
    % Left column: List of modules
    modlist = uicontrol(f, 'Units', 'normalized', ...
                 'Position', [0, 0, lwidth, 1],...
                 'Style', 'listbox',...
                 'String', modules, ...
                 'Callback', @onClickModule);
    % Middle column: Title and description
    infoPanel = uipanel(f, 'Units', 'normalized', 'Position', topPos, ...
                           'TitlePosition', 'centertop', 'fontsize', 18, ...
                           'BorderType', 'none');

    descr = uicontrol(infoPanel, 'Units', 'normalized', ...
                 'Position', [0, .5, 1, .5],...
                 'Style', 'edit',...
                 'min', 0, 'max', 1e3, ...
                 'HorizontalAlignment', 'left');

    expanel = uipanel(infoPanel, 'Units', 'normalized', ...
                      'Position', [0, 0, 1, .5], ...
                      'FontSize', 14, 'TitlePosition', 'centertop',...
                      'BorderType', 'none');
    exinfo = uicontrol(expanel, 'Units', 'normalized', ...
                 'Position', [0, 0, 1, 1],...
                 'Style', 'edit',...
                 'min', 0, 'max', 1e3, ...
                 'HorizontalAlignment', 'left');
    % Reference GUI
    paperPanel = uipanel(f, 'Units', 'normalized', 'Position', endPos);
    paperh = [];
    
    % Right column: List of examples and corresponding buttons
    exlist = uicontrol(f, 'Units', 'normalized', ...
                 'Position', [midwidth + lwidth, .05, lwidth, .95],...
                 'Callback', @setMainBoxText, ...
                 'Style', 'listbox');

    uicontrol(f, 'Units', 'normalized', ...
                 'Position', [midwidth + lwidth, 0, lwidth/3, .05],...
                 'Style', 'pushbutton',...
                 'String', 'Edit', ...
                 'Callback', @onClickExample);
             
    uicontrol(f, 'Units', 'normalized', ...
                 'Position', [midwidth + 4/3*lwidth, 0, lwidth/3, .05],...
                 'Style', 'pushbutton',...
                 'Tag',    'publish', ...
                 'String', 'Publish', ...
                 'Callback', @onClickExample);

     uicontrol(f, 'Units', 'normalized', ...
                 'Position', [midwidth + 5/3*lwidth, 0, lwidth/3, .05],...
                 'Style', 'pushbutton',...
                 'Tag',    'viewhtml', ...
                 'String', 'View HTML', ...
                 'Callback', @onClickExample);
    function onClickModule(src, event)
        ix = get(src, 'Value');
        mods = get(src, 'String');
        module = mods{ix};
        if strcmpi(module, 'mrst')
            % Core mode
            module = 'core';
            name = 'MRST Core';
        else
            name = ['Module "', module, '"'];
        end
        if ishandle(paperh)
            delete(paperh);
        end
        set(infoPanel, 'Title', name);
        
        [examples, examplenames] = getExamples(module);
        set(exlist, 'String', examplenames);
        set(exlist, 'Value', 1);
        paperh = mrstReferencesGUI({module}, paperPanel, [0, 0, 1, 1], 'Relevant literature');
        setMainBoxText();
    end

    function onClickExample(src, event)
        ix = get(exlist, 'Value');
        if isempty(examples)
            return
        end
        fn = examples{ix};
        switch lower(get(src, 'Tag'))
            case 'publish'
                modNo = get(modlist, 'Value');
                if ~strcmpi(modules{modNo}, 'mrst')
                    mrstModule('add', modules{modNo});
                end
                webfile = publish(fn);
                web(webfile);
            case 'viewhtml'
                [fldr, nm] = fileparts(fn);
                expth = fullfile(fldr, 'html', [nm, '.html']);
                if exist(expth, 'file')
                    web(expth);
                else
                    errordlg('Example not published! Run publish first.');
                end
            otherwise
                edit(fn);
        end
    end

    function setMainBoxText(src, event)
        ix = get(exlist, 'Value');
        if isempty(examples)
            fn = '';
        else
            fn = examples{ix};
        end
        
        modNo = get(modlist, 'Value');
        module = modules{modNo};
        
        if isempty(fn)
            exDescr = '';
            sn = '';
        else
            exDescr = help(fn);
            [~, sn] = fileparts(fn);
        end
        modDescr = getDescription(module);
        set(descr,  'String', modDescr);
        set(exinfo, 'String', exDescr);
        set(expanel, 'Title', sn);
    end

    onClickModule(modlist, {'core'});
    
    if nargout > 0
        varargout{1} = f;
    end
    set(f, 'HandleVisibility',  'off')
end

function d = getDescription(module)
    if strcmpi(module, 'mrst')
        pth = ROOTDIR();
    else
        pth = mrstPath(module);
    end
    files = dir(pth);
    d = ['No README available for ', module];

    files = { files(~ [files.isdir]).name };

    i = strcmpi(files, 'readme') | ...
        strcmpi(files, 'readme.txt');

    if any(i),
       readme = fileread(fullfile(pth, files{find(i, 1, 'first')}));

       if ~ isempty(readme),
          d = mrstSplitText(readme);
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

function pos = getPanelPosition()
% Utility to set dimensions of the main panel
    % First, get dimensions of first monitor
    monpos = get(0, 'MonitorPositions');
    monsz = monpos(1, 3:4);
    % Try to create a window that is reasonably large, but not so large
    % that it is bigger than the monitor itself
    w = min(0.8*monsz(1), 1500);
    h = min(0.8*monsz(2), 1000);
    % Place the window centered
    start = (monsz - [w, h])/2;
    pos = [start, w, h];
end
