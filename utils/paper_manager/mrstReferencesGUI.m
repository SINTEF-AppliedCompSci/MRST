function varargout = mrstReferencesGUI(modules, parent, pos)
% Draw a panel for exploring relevant papers for MRST
%
% SYNOPSIS:
%   mrstReferencesGUI();
%   % Embed in an existing user interface
%   mrstReferencesGUI(parent, [0, 0, 1, .25])
%
% DESCRIPTION:
%   Create a panel for exploring publications that use MRST.
%
% PARAMETERS:
%  modules - A cell array containing a list of modules for which
%            publications are to be requested. If "modules" is empty, or
%            not given, it will default to all known publications.
%
%  parent  - Handle to a UI component where reference GUI should be added.
%            If not given, this function will default to a brand new
%            figure.
%
%  pos     - Position argument for embedding in parent. Defaults to 
%            [0, 0, 1, 1], filling the entire surface of the parent.
%
% RETURNS:
%   panel  - Handle to the panel that has been added to the workspace.
%
%   parent - Handle to the parent, or the figure that was created if no
%            parent was given.
% 
% NOTE:
%  Note that all arguments are optional. Passing no input arguments
%
% SEE ALSO:
%   mrstExploreModules, getAvailablePapers

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
    if nargin < 2
        parent = figure;
        set(parent,'Name',          'MRST paper explorer', ...
                   'Toolbar',       'none', ...
                   'NumberTitle',   'off', ...
                   'Units', 'normalized', ...
                   'MenuBar',       'none');
    end
    assert(ishandle(parent), 'Parent handle not valid.');
    
    if nargin < 3
        pos = [0, 0, 1, 1];
    end
    panel = uipanel(parent, 'position', pos);
    papers = getAvailablePapers();
    
    if nargin > 0 && ~isempty(modules)
        if ischar(modules)
            modules = {modules};
        end
        act = arrayfun(@(x) any(ismember(x.modules, modules)), papers);
        papers = papers(act);
    end
    if isempty(papers)
        papers = createPaperStruct('name', 'No papers found');
    end
    
    names = arrayfun(@(x) x.name, papers, 'UniformOutput', false);
    picker = uicontrol(panel, 'Units', 'normalized', ...
                 'Position', [0, .4, 1, .6],...
                 'Style', 'listbox',...
                 'String', names, ...
                 'Callback', @onClickPapers);
             
    txt = uicontrol(panel, 'Units', 'normalized', ...
                 'Position', [0, 0.1, 1, .3],...
                 'Min', 0, 'Max', inf, ...
                 'Style', 'text');
    function onClickPapers(src, event)
        ix = get(src, 'Value');
        
        paper = papers(ix);
        
        if isempty(paper.authors)
            auth = 'No authors';
        else
            auth = ['by ', paper.authors];
        end
        set(txt, 'String', {paper.name, ...
                           auth, ...
                           [paper.published, ' ', formatYear(paper.year)], ...
                           });
    end

    function openFieldBrowser(src, event, fld)
        ix = get(src, 'Value');
        paper = papers(ix);
        f = paper.(fld);
        if isempty(f)
            errordlg('No URL supplied!')
        else
            web(f);
        end
    end
    bg = uipanel(panel, ...
                 'Units', 'normalized', ...
                 'Position', [0, 0, 1, .1]);
    
    uicontrol(bg, 'Units', 'normalized',...
                  'Position', [0, 0, .2, 1], ...
                  'Style', 'pushbutton', ...
                  'String', 'View', ...
                  'callback', @(src, event) openFieldBrowser(src, event, 'url'))
              
    uicontrol(bg, 'Units', 'normalized',...
                  'Position', [.2, 0, .2, 1], ...
                  'Style', 'pushbutton', ...
                  'String', 'Preprint', ...
                  'callback', @(src, event) openFieldBrowser(src, event, 'fileurl'))

    onClickPapers(picker, []);
    % Return panel handle
    if nargout > 0
        varargout{1} = panel;
        if nargout > 1
            varargout{2} = parent;
        end
    end
end

function s = formatYear(yr)
    if yr == -1
        s = '';
    else
        s = ['(', num2str(yr), ')'];
    end
end

