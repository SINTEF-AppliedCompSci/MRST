function mrstExploreModules()
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
    paperh = mrstReferencesGUI(f, midpos);
    
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
        paperh = mrstReferencesGUI(f, midpos, m);
    end

    function onClickExample(src, event)
        ix = get(exlist, 'Value');
        fn = examples{ix};
        if strcmpi(get(src, 'Tag'), 'publishbutton')
            webfile = publish(fn);
            web(webfile);
        else
            edit(fn);
        end
    end

    onClickModule(modlist, []);
end

function d = getDescription(module)
    pth = mrstPath(module);
    files = dir(pth);
    d = 'No description available.';
    for i = 1:numel(files)
        if strcmpi(files(i).name, 'contents.m')
            d = fileread(fullfile(pth, files(i).name));
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

