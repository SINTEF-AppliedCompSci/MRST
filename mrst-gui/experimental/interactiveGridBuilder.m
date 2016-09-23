function varargout = interactiveGridBuilder(varargin)
    opt = struct('image', '');
    opt = merge_options(opt, varargin{:});
    
    if ~isempty(opt.image)
        I = imread(opt.image);
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

    linebutton = uicontrol(panel, 'units', 'normalized', 'Position', [0, 0.8, 1, 0.1],...
                     'String', 'Draw line', 'Style', 'pushbutton', ...
                     'callback', @lineButtonHandler);
   
    pointbutton = uicontrol(panel, 'units', 'normalized', 'Position', [0, 0.7, 1, 0.1],...
                     'String', 'Add point', 'Style', 'pushbutton', ...
                     'callback', @pointButtonHandler);
                 
    undobutton = uicontrol(panel, 'units', 'normalized', 'Position', [0, 0.6, 1, 0.1],...
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
                 
    allbuttons = [drawbutton; linebutton; pointbutton];
    outline = [];
    lines = {[]};
    points = [];
    lastPick = nan;
    
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
        switch lastPick
            case 1
                % Outline
                outline = [];
                msg = 'Grid outline removed!';
            case 2
                % Line
                if numel(lines) > 1
                    lines = lines([1:end-2, end]);
                end
                msg = 'Last line segment removed!';
            case 3
                % Point
                if ~isempty(points)
                    points = points(1:end-1, :);
                end
                msg = 'Last point removed!';
            otherwise
                msg = 'Nothing to remove!';
        end
        displayMessage(msg);
        drawOutline();
    end

    function drawButtonHandler(src, event)
        outline = [];
        set(h, 'WindowButtonDownFcn', @addOutlinePt);
        displayMessage(['Drawing domain outline. Left click to add points'...
            ' and right-click to stop drawing. Polygon is automatically'...
            ' closed on right-click']);
        disableButtons()
    end
    
    function lineButtonHandler(src, event)
        set(h, 'WindowButtonDownFcn', @addLinePt);
        displayMessage(['Adding a line segment. Left click to add points to'...
            ' current line, and right click to mark current line as complete.']);
        disableButtons()
    end

    function pointButtonHandler(src, event)
        set(h, 'WindowButtonDownFcn', @addPoint);
        displayMessage(['Adding points to grid. Left click to add points. '...
            'Right click when done.']);

        disableButtons()
    end

    function drawOutline()
        cla(pax);
        if ~isempty(I)
            imh = imshow(I, 'Parent', pax, 'XData', [0, 1], 'YData', [0, 1]);
            set(imh, 'HitTest', 'off');
            axis normal
        end
        set(gca, 'XAxisLocation', 'bottom');
        set(gca, 'YDir', 'normal');
        hold on
        if ~isempty(outline)
            plot(pax, outline(:, 1), outline(:, 2));
        end
        
        for i = 1:numel(lines)
            if isempty(lines{i})
                continue
            end
            
            plot(lines{i}(:, 1), lines{i}(:, 2), 'r');
        end
        
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
            drawOutline();
            enableButtons();
            return
        end
        
        points = [points; pt];
        lastPick = 3;
        drawOutline();
    end

    function addLinePt(src, event)
        pt = get(h, 'currentpoint');
        pt = getLocalAxisCoordinate(pt);
        
        if ~strcmpi(get(h, 'SelectionType'), 'normal') && ~isempty(lines{end})
            set(h, 'WindowButtonDownFcn', []);
            lines{end+1} = [];
            drawOutline();
            enableButtons();
            return
        end

        if pt(1) > 1 || pt(2) > 1
            return
        end
        
        lines{end} = [lines{end}; pt];
        lastPick = 2;
        drawOutline();
    end

    function addOutlinePt(src, event)
        pt = get(h, 'currentpoint');
        pt = getLocalAxisCoordinate(pt);
        
        if ~strcmpi(get(h, 'SelectionType'), 'normal') && ~isempty(outline)
            outline = [outline; outline(1, :)];
            set(h, 'WindowButtonDownFcn', []);
            drawOutline();
            enableButtons();
            return
        end

        if pt(1) > 1 || pt(2) > 1
            return
        end
        
        outline = [outline; pt];
        lastPick = 1;
        drawOutline();
    end

    function out = packOutput(src, event)
        out = struct();
        out.points = points;
        out.lines = lines(1:end-1);
        if ~isempty(outline)
            out.outline = outline(1:end-1, :);
        else
            out.outline = [];
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


    drawOutline();
    if nargout > 0
        uiwait(gcf);
    end
end

