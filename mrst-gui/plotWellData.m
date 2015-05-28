function [varargout] = plotWellData(G, W, varargin)
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
                 'EdgeColor', {'none'}, ...
                 'linePlot', false);
    opt = merge_options(opt, varargin{:});
    
    if isempty(opt.TextColor)
        opt.TextColor = opt.color;
    end
    nW    = numel(W);
    htop  = zeros([1, nW]);
    htext = zeros([1, nW]);
    hs    = zeros([1, nW]);
    hline = zeros([1, nW]);
    
    dohold = ishold();
    
    hold on
    for i = 1:numel(W)
        if ~opt.linePlot
            [hs(i), htop(i)] = plotSingleWell(G, W(i), data{i}, opt);
        else
            [hs(i), htop(i)] = plotSingleWellTraj(G, W(i), data{i}, opt);
        end
        [htext(i), hline(i)] = plotSingleWellLabel(G, W(i), opt);
    end
    
    
    if ~dohold
        hold off
    end
    if nargout > 0
        varargout{1} = htop;
    end
    
    if nargout > 1
        varargout{2} = htext;
    end
    
    if nargout > 2
        varargout{3} = hs;
    end
    
    if nargout > 3
        varargout{4} = hline;
    end
end

function [hs, htop] = plotSingleWell(G, W, data, opt)
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
    [hs, htop] = plotCurveTube(curve, opt.radius*R, 'refine', 5, 'data', data,...
                                'EdgeColor', opt.EdgeColor,...
                                'color', opt.color, ...
                                'interpMethod', opt.interpMethod);
end

function [hs, htop] = plotSingleWellTraj(G, W, data, opt)%#ok
    if ~isfield(W, 'trajectory')
        traj = G.cells.centroids(W.cells,:);
        warning('Well-trajectory for plotting not found, using connection cells');
    else
        traj = W.trajectory;
    end
    hs   = plot3(traj(:,1), traj(:,2), traj(:,3), 'LineWidth', 4, 'color', opt.color);
    htop = plot3(traj([1,end],1), traj([1,end],2), traj([1,end],3), 'o', 'MarkerFaceColor', opt.color, 'color', opt.color);
end


function [htext, hline] = plotSingleWellLabel(G, W, opt)
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
    
    htext = text(topBore(1), topBore(2), height2, W.name, ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'center', ...
        'Rotation', opt.labelRotation, ...
        'FontSize', opt.labelFontSize,...
        'EdgeColor', ec, ...
        'Interpreter', 'tex', ...
        'Color', opt.TextColor, ...
        'backgroundcolor', opt.labelBackgroundColor);
    hline = plot3(d(:, 1), d(:, 2), d(:, 3), 'LineWidth', 2, 'Color', opt.color);
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
