function panel = mrstReferencesGUI(parent, pos, modules)
    if nargin == 0
        parent = figure;
        set(parent,'Name',          'MRST paper explorer', ...
                   'Toolbar',       'none', ...
                   'NumberTitle',   'off', ...
                   'Units', 'normalized', ...
                   'MenuBar',       'none');
        pos = [0, 0, 1, 1];
    end
    panel = uipanel(parent, 'position', pos);
    
    papers = getAvailablePapers();
    
    if nargin == 3 && ~isempty(modules)
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
end

function s = formatYear(yr)
    if yr == -1
        s = '';
    else
        s = ['(', num2str(yr), ')'];
    end
end

