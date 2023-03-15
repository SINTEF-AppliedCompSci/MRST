function W = addHybridCellsToWell(W, G, rock, varargin)

    opt = struct('aperture', []);
    opt = merge_options(opt, varargin{:});

    % Topology (+1 to avoid zero index)
    N = G.faces.neighbors + 1;
    % Hybrid cell mask (pad with leading false to comply above trick)
    hc  = G.cells.hybrid > 0;
    hc  = [false; hc];
    % Number of matrix cells
    ncm = G.cells.num - nnz(G.cells.hybrid);
    for i = 1:numel(W)
        % Well cell mask
        wc = false(G.cells.num+1, 1); wc(W(i).cells+1) = true;
        % Find well cells that are adjecent to hybrid cells
        ix = any(hc(N), 2) & any(wc(N), 2);
        if any(ix)
            cells  = reshape(N(ix,:), [], 1)-1;
            hcells = unique(cells(cells > ncm));
            W(i)   = addWellCells(W(i), G, rock, hcells, opt.aperture);
        end
    end
    
end

%-------------------------------------------------------------------------%
function W = addWellCells(W, G, rock, cells, aperture)
% Add cells

    nc0     = numel(W.cells);
    W.cells = [W.cells; cells];
    nc      = numel(W.cells);
    % Safeguard against negative values by enlarging cellDims if necessary
    [dx, dy, dz] = cellDims(G, cells);
    if ~isempty(aperture)
        dz = aperture;
        if numel(dz) == G.cells.num, dz = dz(cells); end
    end
    re   = 2*0.14*sqrt(dx.^2 + dy.^2)/2;
    C    = max(1.1*W.r(1)./re,1);
    dx   =  [dx.*C, dy.*C, dz];
    % Compute well index (should probably be calculated in a different way)
    WI   = computeWellIndex(G, rock, W.r(1), cells, 'cellDims', dx);
%     WI   = computeWellIndex(G, rock, W.r(1), cells);
%     WI   = min(WI, max(W.WI*1000));
    W.WI = [W.WI; WI];
    % Add delta z field
    dZ   = getDeltaZ(G, cells, W.refDepth);
    W.dZ = [W.dZ; dZ];
    [~, ix] = sort(W.dZ);
    
    W.cells = W.cells(ix);
    W.WI = W.WI(ix);
    W.dZ = W.dZ(ix);
    
    % Add remaining fields by repeating W.(fn)(1)
    for fn = fieldnames(W)'
        v = W.(fn{1});
        if size(v,1) ~= nc0, continue; end
        W.(fn{1}) = repmat(v(1,:), nc, 1);
    end
    
end

%-------------------------------------------------------------------------%
function dZ = getDeltaZ(G, cells, refDepth)
% Compute distance from ref height to perforation (Copied from addWell)

    direction = gravity();
    dims      = G.griddim;
    if norm(direction(1:dims)) > 0
       direction = direction ./ norm(direction(1:dims));
    else
       direction = zeros(1, dims);
       if dims > 2
          direction(end) = 1;
       end
    end
    xyz = G.cells.centroids(cells, :);
    xyz(isnan(xyz)) = 0;
    Z = xyz * direction(1:dims).';
    dZ = Z - refDepth;
    
end