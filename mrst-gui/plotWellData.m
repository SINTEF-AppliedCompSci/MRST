function plotWellData(G, W, varargin)
    if mod(numel(varargin), 2) == 1
        data = varargin{1};
        varargin = varargin(2:end);
        if ~iscell(data) && numel(W) == 1
            data = {data};
        end
    else
        data = cell(numel(W), 1);
    end
    opt = struct('aspectScale', true, ...
                 'showInactive', false, ...
                 'radialData', false, ...
                 'colorByData', true, ...
                 'labelBackgroundColor', 'none', ...
                 'labelFontSize',   16, ...
                 'labelHeight', .1, ...
                 'labelRotation', 0, ...
                 'radius', 1, ...
                 'interpMethod', 'pchip', ...
                 'color',   {[1 0 0]}, ...
                 'TextColor', [], ...
                 'EdgeColor', {'none'});
    opt = merge_options(opt, varargin{:});
    
    if isempty(opt.TextColor)
        opt.TextColor = opt.color;
    end
    
    dohold = ishold();
    
    hold on
    for i = 1:numel(W)
        plotSingleWell(G, W(i), data{i}, opt);
        plotSingleWellLabel(G, W(i), opt);
    end
    
    if ~dohold
        hold off
    end
end

function plotSingleWell(G, W, data, opt)
    [R, curve] = computeEffectiveRadius(G, W, opt);
    
    if opt.radialData && numel(data) > 1
        d = data - min(data);
        d = d./max(d);
        
        R = repmat(2*(d*0.9 + 0.1), 1, size(R, 2)).*R;
    end
    
    if isfield(W, 'cstatus') && ~opt.showInactive
        disabled = ~W.cstatus;
        R(disabled, :) = repmat(0.001*min(R), sum(disabled > 0), 1);
    end
    
    if ~opt.colorByData
        data = [];
    end
    plotCurveTube(curve, opt.radius*R, 'refine', 5, 'data', data,...
                                'EdgeColor', opt.EdgeColor,...
                                'color', opt.color, ...
                                'interpMethod', opt.interpMethod);
end

function plotSingleWellLabel(G, W, opt)
    if G.griddim == 3
        top = min(G.nodes.coords(:, 3));
        dist = max(G.nodes.coords(:, 3)) - top;
        
        topBore = G.cells.centroids(W.cells(1), :);
    else
        top = 0;
        dist = 1;
        topBore = [G.cells.centroids(W.cells(1), :), 0];
    end
    height1 = top - 0.95*opt.labelHeight*dist;
    height2 = top - opt.labelHeight*dist;
    
    d = [topBore; topBore];
    d(2, 3) = height1;
    
    if strcmpi(opt.labelBackgroundColor, 'none')
        ec = 'none';
    else
        ec = opt.color;
    end
    
    text(topBore(1), topBore(2), height2, W.name, ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'center', ...
        'Rotation', opt.labelRotation, ...
        'FontSize', opt.labelFontSize,...
        'EdgeColor', ec, ...
        'Interpreter', 'tex', ...
        'Color', opt.TextColor, ...
        'backgroundcolor', opt.labelBackgroundColor)
    plot3(d(:, 1), d(:, 2), d(:, 3), 'LineWidth', 2, 'Color', opt.color);
end

function [R, curve, cells] = computeEffectiveRadius(G, W, opt)
    cells = W.cells;
    curve = G.cells.centroids(cells, :);
    if opt.aspectScale
        aspect =  max(G.nodes.coords) - min(G.nodes.coords);
        aspect = aspect./sum(aspect);
    else
        aspect = ones(1, G.griddim);
    end
    r = ones(numel(cells), 1);
    
    % Magic scaling factor
    sfac = sqrt(mean(G.faces.areas));
    R = bsxfun(@mtimes, r, aspect)*sfac;
end
