function G = convertBuilderToPEBI(out, n, varargin)
    require upr
    fault = out.lines;
    bdr = out.outline;
    wells = cell(1, size(out.points, 1));
    for i = 1:numel(wells)
        wells{i} = out.points(i, :);
    end
    
    L = fault;
    % Trick the gridder to make outline of the domain as faces so we can
    % remove the final cells
    if ~isempty(bdr)
        bdr_closed = [bdr; bdr(1, :)];
        for i = 1:size(bdr, 1);
            L = [L, bdr_closed([i, i+1], :)];
        end
    end
    G = pebiGrid(1/sqrt(n), [1, 1], ...
        'faultLines', L, 'wellLines', wells, ...
        varargin{:});
    
    if ~isempty(bdr)
        % Extract subgrid defined by bounding box
        G = computeGeometry(G);
        keep = inpolygon(G.cells.centroids(:, 1), G.cells.centroids(:, 2),...
                         bdr_closed(:, 1), bdr_closed(:, 2));
        G = extractSubgrid(G, keep);
    end
end
