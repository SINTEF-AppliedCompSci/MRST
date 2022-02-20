function h = plotWellData(G, W, varargin)
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
                 'labelBackgroundColor', {{'none'}}, ...
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
    nW     = numel(W);
    h = repmat(struct('label', [], 'connector', [], 'body', []), [1, nW]);
    dohold = ishold();
    hold on
    for i = 1:nW
        if ~opt.linePlot
            if ~strcmp(W(i).type, 'aquifer')
                [h(i).body] = plotSingleWell(G, W(i), data{i}, opt); 
            end
        else
            h(i).body = plotSingleWellTraj(G, W(i), data{i}, opt);
        end
        [h(i).label, h(i).connector] = plotSingleWellLabel(G, W(i), opt);
    end
    
    
    if ~dohold
        hold off
    end
end

function body = plotSingleWell(G, W, data, opt)
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
    body = [hs, htop];
end

function body = plotSingleWellTraj(G, W, data, opt)%#ok
    if ~isfield(W, 'trajectory')
        traj = G.cells.centroids(W.cells,:);
        warning('Well-trajectory for plotting not found, using connection cells');
    else
        traj = W.trajectory;
    end
    body  = line(traj(:,1), traj(:,2), traj(:,3), 'LineWidth', 5, 'color', opt.color);
    %htop = plot3(traj([1,end],1), traj([1,end],2), traj([1,end],3), 'o', 'MarkerFaceColor', opt.color, 'color', opt.color);
end


function [htext, hline] = plotSingleWellLabel(G, W, opt)
    if G.griddim == 3
        if isfield(G,'nodes')
           top = min(G.nodes.coords(:, 3));
           dist = max(G.nodes.coords(:, 3)) - top;
        else
           top = min(G.cells.centroids(:, 3));
           dist = max(G.cells.centroids(:, 3)) - top;
        end
        topBore = G.cells.centroids(W.cells(1), :);
    else
        top = 0;
        dist = 1;
        topBore = [G.cells.centroids(W.cells(1), :), 0];
    end
    height1 = top - 0.9*opt.labelHeight*dist;
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
        'backgroundcolor', opt.labelBackgroundColor{1});
    hline = plot3(d(:, 1), d(:, 2), d(:, 3), 'LineWidth', 2, 'Color', opt.color);
end

function [R, curve, cells] = computeEffectiveRadius(G, W, opt)
    cells = W.cells;
    curve = G.cells.centroids(cells, :);
    if opt.aspectScale && isfield(G,'nodes')
        aspect =  max(G.nodes.coords) - min(G.nodes.coords);
        aspect = aspect./sum(aspect);
    elseif isfield(G,'parent')
        aspect =  max(G.parent.nodes.coords) - min(G.parent.nodes.coords);
        aspect = aspect./sum(aspect);
    else
        aspect = ones(1, G.griddim);
    end
    r = ones(numel(cells), 1);
    
    % Magic scaling factor
    sfac = sqrt(mean(G.faces.areas));
    R = bsxfun(@mtimes, r, aspect)*sfac;
end
