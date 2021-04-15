function varargout = mrstReferencesGUI(modules, parent, pos, name)
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
%   `mrstExploreModules`, `getAvailablePapers`

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
    if nargin < 4
        name = '';
    end
    panel = uipanel(parent, 'position', pos, 'Title', name, ...
                            'BorderType', 'none', ...
                            'TitlePosition', 'centertop', ...
                            'FontSize', 14);
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
                 'Min', 0, 'Max', 1e3, ...
                 'HorizontalAlignment', 'left',...
                 'Style', 'text');
    function onClickPapers(src, event)
        ix = get(src, 'Value');
        paper = papers(ix);
        
        if isempty(paper.authors)
            auth = ' No authors';
        else
            auth = [' by ', paper.authors];
        end
        
        doi = paper.doi;
        if ~isempty(doi)
            doi = [' DOI: ', doi];
        end
        set(txt, 'String', [paper.name, ...
                           auth, ...
                           [' ', paper.published, ' ', formatYear(paper.year)], ...
                           doi,...
                           ]);
        if isempty(paper.url)
            set(openurl, 'Enable', 'off');
        else
            set(openurl, 'Enable', 'on');
        end
        if isempty(paper.fileurl)
            set(openfile, 'Enable', 'off');
        else
            set(openfile, 'Enable', 'on');
        end
        if isempty(paper.doi)
            set(getcitation, 'Enable', 'off');
        else
            set(getcitation, 'Enable', 'on');
        end
    end

    function openFieldBrowser(src, event, fld)
        ix = get(picker, 'Value');
        paper = papers(ix);
        f = paper.(fld);
        if isempty(f)
            errordlg('No URL supplied!')
        else
            web(f, '-browser');
        end
    end

    function getCitation(src, event)
        ix = get(picker, 'Value');
        paper = papers(ix);
        f = paper.doi;
        if isempty(f)
            errordlg('No DOI supplied for paper!')
        else
            % We query the DOI database and request bibtex as output
            if exist('weboptions', 'file')
                opt = weboptions('ContentType', 'text', ...
                                 'KeyName', 'Accept', ...
                                 'KeyValue', 'application/x-bibtex');
                bib = webread(['http://doi.org/', f], opt);
                bib = regexprep(bib,' +',' ');
                ch = figure('Position',      getCitationPosition(), ...
                            'Name',          paper.name, ...
                            'Toolbar',       'none', ...
                            'NumberTitle',   'off', ...
                            'MenuBar',       'none');
                wrn = ['BibTeX entry automatically generated from DOI.org. ',...
                       'Please verify information before using!'];
                uicontrol(ch, 'Style', 'edit', ...
                              'String', 'bibtex', ...
                              'Max', 1e3, ...
                              'String', {wrn, '', bib}, ...
                              'Units', 'normalized', ...
                              'HorizontalAlignment', 'left', ...
                              'position', [0, 0, 1, 1]);
            else
                % Fallback to DOI2BIB website
                web(['http://www.doi2bib.org/#/doi/', f]);
            end
        end
    end

    bg = uipanel(panel, ...
                 'Units', 'normalized', ...
                 'Position', [0, 0, 1, .1]);
    
    openurl = uicontrol(bg, 'Units', 'normalized',...
                  'Position', [0, 0, .2, 1], ...
                  'Style', 'pushbutton', ...
                  'String', 'View', ...
                  'callback', @(src, event) openFieldBrowser(src, event, 'url'));
              
    openfile = uicontrol(bg, 'Units', 'normalized',...
                  'Position', [.2, 0, .2, 1], ...
                  'Style', 'pushbutton', ...
                  'String', 'Preprint', ...
                  'callback', @(src, event) openFieldBrowser(src, event, 'fileurl'));

    getcitation = uicontrol(bg, 'Units', 'normalized',...
                  'Position', [.4, 0, .2, 1], ...
                  'Style', 'pushbutton', ...
                  'String', 'Export citation', ...
                  'callback', @(src, event) getCitation(src, event));

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

function pos = getCitationPosition()
    % First, get dimensions of first monitor
    monpos = get(0, 'MonitorPositions');
    monsz = monpos(1, 3:4);
    % Try to create a window that is reasonably large, but not so large
    % that it is bigger than the monitor itself
    w = min(0.5*monsz(1), 1000);
    h = min(0.5*monsz(2), 400);
    % Place the window centered
    start = (monsz - [w, h])/2;
    pos = [start, w, h];
end
