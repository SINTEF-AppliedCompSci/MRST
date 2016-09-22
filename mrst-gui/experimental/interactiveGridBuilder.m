function varargout = interactiveGridBuilder()
    h = figure();
    
    axwidth = 0.75;
    pax = subplot('Position', [0, 0, axwidth, 1]);
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
    if nargout > 0
        returnbutton = uicontrol(panel, 'units', 'normalized', 'Position', [0, 0.0, 1, 0.1],...
                         'String', 'Return values', 'Style', 'pushbutton', ...
                         'callback', @returnButtonHandler);
    end
    allbuttons = [drawbutton; linebutton; pointbutton];
    outline = [];
    lines = {[]};
    points = [];
    
    function disableButtons()
        set(allbuttons, 'Enable', 'off');
    end

    function enableButtons()
        set(allbuttons, 'Enable', 'on');
    end
    function drawButtonHandler(src, event)
        outline = [];
        set(h, 'WindowButtonDownFcn', @addOutlinePt);
        disableButtons()
    end
    
    function lineButtonHandler(src, event)
        set(h, 'WindowButtonDownFcn', @addLinePt);
        disableButtons()
    end

    function pointButtonHandler(src, event)
        set(h, 'WindowButtonDownFcn', @addPoint);
        disableButtons()
    end

    function drawOutline()
        cla(pax);
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
        figpos = get(h, 'Position');

        pt = pt./figpos(3:end);
        pt(1) = pt(1)/axwidth;
        
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
        drawOutline();
    end

    function addLinePt(src, event)
        pt = get(h, 'currentpoint');
        figpos = get(h, 'Position');
        
        if ~strcmpi(get(h, 'SelectionType'), 'normal') && ~isempty(lines{end})
            set(h, 'WindowButtonDownFcn', []);
            lines{end+1} = [];
            drawOutline();
            enableButtons();
            return
        end
        pt = pt./figpos(3:end);
        pt(1) = pt(1)/axwidth;
        
        if pt(1) > 1 || pt(2) > 1
            return
        end
        
        lines{end} = [lines{end}; pt];
        drawOutline();
    end

    function addOutlinePt(src, event)
        pt = get(h, 'currentpoint');
        figpos = get(h, 'Position');
        
        if ~strcmpi(get(h, 'SelectionType'), 'normal') && ~isempty(outline)
            outline = [outline; outline(1, :)];
            set(h, 'WindowButtonDownFcn', []);
            drawOutline();
            enableButtons();
            return
        end
        pt = pt./figpos(3:end);
        pt(1) = pt(1)/axwidth;
        
        if pt(1) > 1 || pt(2) > 1
            return
        end
        
        outline = [outline; pt];
        drawOutline();
    end

    function returnButtonHandler(src, event)
        out = struct();
        out.points = points;
        out.lines = lines(1:end-1);
        if ~isempty(outline)
            out.outline = outline(1:end-1, :);
        else
            out.outline = [];
        end
        
        varargout{1} = out;
        uiresume();
    end
    if nargout > 0
        uiwait(gcf);
    end
end