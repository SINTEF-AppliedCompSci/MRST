function varargout = interactiveGridBuilder(varargin)

    if mod(nargin, 2) ~= 0
        input = varargin{1};
        varargin = varargin(2:end);
    else
        input = [];
    end
    opt = struct('image', '');
    opt = merge_options(opt, varargin{:});
    
    if ~isempty(opt.image)
        if isnumeric(opt.image)
            I = opt.image;
        else
            I = imread(opt.image);
        end
        I = flipud(I);
    else
        I = [];
    end
    h = figure();
    
    axwidth = 0.75;
    axheight = 0.8;
    pax = subplot('Position', [0, 1-axheight, axwidth, axheight]);
    panel = uipanel(h, 'Position', [axwidth, 0, 1-axwidth, 1]);
    set(pax,'HitTest','off')
    

    
    drawbutton = uicontrol(panel, 'units', 'normalized', 'Position', [0, 0.9, 1, 0.1],...
                     'String', 'Draw outline', 'Style', 'pushbutton', ...
                     'callback', @drawButtonHandler);

    faultbutton = uicontrol(panel, 'units', 'normalized', 'Position', [0, 0.8, 1, 0.1],...
                     'String', 'Draw fault', 'Style', 'pushbutton', ...
                     'callback', @(src, event) lineButtonHandler(src, event, 'fault'));

    wellbutton = uicontrol(panel, 'units', 'normalized', 'Position', [0, 0.7, 1, 0.1],...
                     'String', 'Draw well', 'Style', 'pushbutton', ...
                     'callback', @(src, event) lineButtonHandler(src, event, 'well'));

    pointbutton = uicontrol(panel, 'units', 'normalized', 'Position', [0, 0.6, 1, 0.1],...
                     'String', 'Add point', 'Style', 'pushbutton', ...
                     'callback', @pointButtonHandler);
                 
    undobutton = uicontrol(panel, 'units', 'normalized', 'Position', [0, 0.5, 1, 0.1],...
                     'String', 'Undo', 'Style', 'pushbutton', ...
                     'callback', @undoButtonHandler);
    if nargout > 0
        returned = false;
        returnbutton = uicontrol(panel, 'units', 'normalized', 'Position', [0, 0.0, 1, 0.1],...
                         'String', 'Return values', 'Style', 'pushbutton', ...
                         'callback', @(src, event) returnButtonHandler(src, event, false));
        set(gcf, 'CloseRequestFcn', @(varargin) returnButtonHandler(h, [], true));
        offset = 0.1;
    else
        offset = 0;
    end
    uicontrol(panel, 'units', 'normalized', 'Position', [0, offset, 1, 0.1],...
                     'String', 'Save in workspace', 'Style', 'pushbutton', ...
                     'callback', @workspaceButtonHandler);
                 
    msgbox = uicontrol(h, 'units', 'normalized', 'Position', [0, 0, axwidth, 1 - axheight],...
                     'String', '', 'Style', 'edit', 'Min', 0, 'Max', 1e6);
                 
    allbuttons = [drawbutton; faultbutton; wellbutton; pointbutton; undobutton];
    outline = [];
    faults = {[]};
    wells = {[]};
    points = [];
    lastPick = [];
    unpackOutput(input);
    
    function displayMessage(msg)
        set(msgbox, 'String', msg);
    end
    
    function disableButtons()
        set(allbuttons, 'Enable', 'off');
    end

    function enableButtons()
        set(allbuttons, 'Enable', 'on');
        displayMessage('');
    end

    function undoButtonHandler(src, event)
        if isempty(lastPick)
            msg = 'Nothing to remove';
        else
            switch lastPick(end)
                case 1
                    % Outline
                    outline = [];
                    msg = 'Grid outline removed!';
                case 2
                    % Line
                    if numel(faults) > 1
                        faults = faults([1:end-2, end]);
                    end
                    msg = 'Last fault removed!';
                case 3
                    % Line
                    if numel(wells) > 1
                        wells = wells([1:end-2, end]);
                    end
                    msg = 'Last well trajectory removed!';
                case 4
                    % Point
                    if ~isempty(points)
                        points = points(1:end-1, :);
                    end
                    msg = 'Last point removed!';
                otherwise
                    error('Internal error: Unknown ID for undo');
            end
            lastPick = lastPick(1:end-1);
        end
        displayMessage(msg);
        drawAllElementsOnAxis();
    end

    function drawButtonHandler(src, event)
        outline = [];
        set(h, 'WindowButtonDownFcn', @addOutlinePt);
        displayMessage(['Drawing domain outline. Left click to add points'...
            ' and right-click to stop drawing. Polygon is automatically'...
            ' closed on right-click']);
        disableButtons()
    end
    
    function lineButtonHandler(src, event, type)
        set(h, 'WindowButtonDownFcn', @(src, event) addLinePt(src, event, type));
        displayMessage(['Adding a ', type, ' segment. Left click to add points to'...
            ' current line, and right click to mark current line as complete.']);
        disableButtons()
    end

    function pointButtonHandler(src, event)
        set(h, 'WindowButtonDownFcn', @addPoint);
        displayMessage(['Adding points to grid. Left click to add points. '...
            'Right click when done.']);

        disableButtons()
    end

    function drawAllElementsOnAxis()
        cla(pax);
        % Draw background image, if it exists
        if ~isempty(I)
            imh = imshow(I, 'Parent', pax, 'XData', [0, 1], 'YData', [0, 1]);
            set(imh, 'HitTest', 'off');
            axis normal
        end
        set(gca, 'XAxisLocation', 'bottom');
        set(gca, 'YDir', 'normal');
        hold on
        % Draw outline
        if ~isempty(outline)
            plot(pax, outline(:, 1), outline(:, 2), 'k', 'linewidth', 2);
        end
        % Draw the faults
        for i = 1:numel(faults)
            if isempty(faults{i})
                continue
            end
            plot(faults{i}(:, 1), faults{i}(:, 2), 'r', 'linewidth', 2);
        end
        % Draw well trajectories
        for i = 1:numel(wells)
            if isempty(wells{i})
                continue
            end
            plot(wells{i}(:, 1), wells{i}(:, 2), 'b', 'linewidth', 2);
        end
        % Draw the points
        if ~isempty(points)
            plot(points(:,1), points(:, 2), 'bo');
        end

        ylim([0, 1]);
        xlim([0, 1]);
    end

    function addPoint(src, event)
        pt = get(h, 'currentpoint');
        pt = getLocalAxisCoordinate(pt);
        
        if pt(1) > 1 || pt(2) > 1
            return
        end
        
        if ~strcmpi(get(h, 'SelectionType'), 'normal')
            set(h, 'WindowButtonDownFcn', []);
            drawAllElementsOnAxis();
            enableButtons();
            return
        end
        
        points = [points; pt];
        lastPick = [lastPick; 4];
        drawAllElementsOnAxis();
    end

    function addLinePt(src, event, type)
        pt = get(h, 'currentpoint');
        pt = getLocalAxisCoordinate(pt);
        isFault = strcmpi(type, 'fault');
        if ~strcmpi(get(h, 'SelectionType'), 'normal') && ...
            ((~isempty(faults{end}) && isFault) || (~isempty(wells{end}) && ~isFault))
            set(h, 'WindowButtonDownFcn', []);
            
            if isFault
                if ~isempty(faults{end})
                    faults{end+1} = [];
                    lastPick = [lastPick; 2];
                end
            else
                if ~isempty(wells{end})
                    wells{end+1} = [];
                    lastPick = [lastPick; 3];
                end
            end
            drawAllElementsOnAxis();
            enableButtons();
            return
        end

        if pt(1) > 1 || pt(2) > 1
            return
        end
        if isFault
            faults{end} = [faults{end}; pt];
        else
            wells{end} = [wells{end}; pt];
        end
        
        drawAllElementsOnAxis();
    end

    function addOutlinePt(src, event)
        pt = get(h, 'currentpoint');
        pt = getLocalAxisCoordinate(pt);
        
        if ~strcmpi(get(h, 'SelectionType'), 'normal') && ~isempty(outline)
            outline = [outline; outline(1, :)];
            set(h, 'WindowButtonDownFcn', []);
            drawAllElementsOnAxis();
            enableButtons();
            lastPick = [lastPick; 1];
            return
        end

        if pt(1) > 1 || pt(2) > 1
            return
        end
        
        outline = [outline; pt];
        drawAllElementsOnAxis();
    end

    function out = packOutput(src, event)
        out = struct();
        out.points = points;
        out.faults = faults(1:end-1);
        out.wells = wells(1:end-1);
        if ~isempty(outline)
            out.outline = outline(1:end-1, :);
        else
            out.outline = [];
        end
    end

    function unpackOutput(input)
        if isempty(input)
            return
        end
        points = input.points;
        faults = [input.faults, {[]}];
        wells = [input.wells, {[]}];
        points = input.points;
        outline = input.outline;
        if ~isempty(outline)
            outline = [outline; outline(1, :)];
        end
    end

    function workspaceButtonHandler(src, event)
        out = packOutput();
        nm = 'GridBuilderOutput';
        assignin('base', nm, out);
        displayMessage(['Stored output as variable ''', nm, ...
                        ''' in workspace.']);
    end

    function returnButtonHandler(src, event, wantClose)
        if returned
            if wantClose
                delete(h)
            end
            return;
        end
        varargout{1} = packOutput();
        returned = true;
        uiresume();
        set(returnbutton, 'Enable', 'off');
        if wantClose
            delete(h)
        end
    end
    function pt = getLocalAxisCoordinate(pt)
        figpos = get(h, 'Position');
        pt = pt./figpos(3:end);
        pt(2) = pt(2) - (1-axheight);
        pt = pt./[axwidth, axheight];
    end


    drawAllElementsOnAxis();
    if nargout > 0
        uiwait(gcf);
    end
end

