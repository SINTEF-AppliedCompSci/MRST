function varargout = gridSliceDice(G, varargin)
hasOutput = nargout > 0;

opt = struct('Location', [0 0 1 1],...
             'Parent',   [], ...
             'Callback', [], ...
             'getHandle', []);

opt = merge_options(opt, varargin{:});

if isempty(opt.Parent)
    parent = figure('Toolbar','none', 'MenuBar', 'none');
    ph = uipanel(parent, 'Position', [0 .8 1 .2]);
else
    parent = opt.Parent;
    ph = uipanel(parent, 'Position', opt.Location);
end

    columnname =   {'Type', 'Height', 'Width', 'Radius', 'X', 'Y', 'Active'};
    columnformat = {{'Square', 'Slice', 'Rectangle'}, 'numeric', 'numeric', 'numeric', 'numeric', {'On', 'Off'}};
    columneditable =  [false true true true];

    sliceh =  uitable('Units','normalized',...
                      'Position', [0 0 1 .8], ...
                      'Data', {}, ...
                      'ColumnName', columnname,...
                      'ColumnFormat', columnformat,...
                      'ColumnEditable', columneditable, ...
                      'CellEditCallback', @applySelection);

    uicontrol(ph, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'String', 'Square', ...
                   'Position', [0 0 .2 1],...
                   'Callback', @addHandler);

    uicontrol(ph, 'Style', 'pushbutton',...
                   'Units', 'normalized',...
                   'String', 'Circle', ...
                   'Position', [.2 0 .2 1],...
                   'Callback', @addHandler);



function addHandler(src, event)
    seltype = get(src, 'String');
    callback = @(start, stop) selectCallback(start, stop, seltype);
    h = opt.getHandle();
    view(get(h, 'Parent'), 0, 90)
    set(h, 'ButtonDownFcn', @(src, event) interactiveSelection(seltype, callback, opt.getHandle))
end

function applySelection(src, event)
    selection = getSelection();
    if hasOutput
        varargout{1} = selection;
        close(parent)
    end
    if ~isempty(opt.Callback)
        opt.Callback(selection)
    end
end

function selectCallback(start, stop, type)
    Data = get(sliceh, 'Data');
    start = start(1,1:2);
    mid = (start + stop)/2;

    switch lower(type)
        case 'square'
            h = abs(stop(1) - start(1));
            w = abs(stop(2) - start(2));

            Data = [Data; {type, h, w, NaN, mid(1), mid(2), 'On'}];
        case 'circle'
            r = norm(stop - start)/2;
            Data = [Data; {type, NaN, NaN, r, mid(1), mid(2), 'On'}];
    end
    set(sliceh, 'Data', Data)
    applySelection(sliceh, []);
end

function selection = getSelection()
    Data = get(sliceh, 'Data');
    x = G.cells.centroids(:,1);
    y = G.cells.centroids(:,2);

    selection = true(G.cells.num, 1);
    for i = 1:size(Data, 1);
        switch lower(Data{i, 1})
            case 'square'
                selection = selection & ...
                            abs(x - Data{i, 5}) <= Data{i, 2}/2 & ...
                            abs(y - Data{i, 6}) <= Data{i, 3}/2;
            case 'circle'
                selection = selection & ...
                            (x - Data{i, 5}).^2 + (y - Data{i, 6}).^2 < Data{i, 4}^2;
        end
    end
end

end
