function G = assembleGrids(Gs)
% Assemble multiple grids, but does not merge common faces and does not 
% handle boundary intersections
%
% SYNOPSIS:
%   Gf = assembleGrids(Gs)
%
% PARAMETERS:
%   Gs  - Grids, nGrid x 1, cell
%
% RETURNS:
%   G - The combined grid
%
% EXAMPLE:
%  G1 = cartGrid([20, 10], [20, 10]);
%  G1 = computeGeometry(G1);
%  cC1 = find(G1.cells.centroids(:,1) < 6 & G1.cells.centroids(:,1) > 3 & ...
%      G1.cells.centroids(:,2) < 6 & G1.cells.centroids(:,2) > 3);
%  cC2 = find(G1.cells.centroids(:,1) < 18 & G1.cells.centroids(:,1) > 14 & ...
%      G1.cells.centroids(:,2) < 4 & G1.cells.centroids(:,2) > 2);
%  G1 = removeCells(G1, [cC1; cC2]);
% 
%  G2 = cartGrid([9, 9], [3, 3]);
%  G2.nodes.coords(:,1) = G2.nodes.coords(:,1) + 3;
%  G2.nodes.coords(:,2) = G2.nodes.coords(:,2) + 3;
% 
%  G3 = cartGrid([12, 6], [4, 2]);
%  G3.nodes.coords(:,1) = G3.nodes.coords(:,1) + 14;
%  G3.nodes.coords(:,2) = G3.nodes.coords(:,2) + 2;
% 
%  figure, hold on; axis equal off
%  plotGrid(G1, 'facecolor', [0, 113, 188]/255)
%  plotGrid(G2, 'facecolor', [216, 82, 24]/255)
%  plotGrid(G3, 'facecolor', [118, 255, 47]/255)
%
% SEE ALSO:
%   `glue2DGrid`

    Gs = convertToColumn(Gs);
    assert(size(Gs, 2) == 1, 'Grid cellarray must be one dimensional!')
    
    % Create an empty struct
    G = struct;

    % G.cells -------------------------------------------------------------
    % num --------------------------------
    nc = cellfun(@(x)x.cells.num, Gs);
    G.cells.num = sum(nc);

    % facePos ----------------------------
    ncf = cellfunUniOut(@(x)diff(x.cells.facePos), Gs);
    ncf = cell2mat(ncf);
    G.cells.facePos = cumsum([1; ncf]);

    % faces ------------------------------
    nf  = cellfun(@(x)x.faces.num, Gs);
    nf_cumsum = cumsum([0; nf]);
    nf_cumsum = mat2cell(nf_cumsum, ones(length(nf_cumsum), 1));
    cf  = cellfunUniOut(@(x)x.cells.faces(:,1), Gs);
    cf  = cellfun(@(x, y) x + y, cf, nf_cumsum(1:end-1), 'UniformOutput', false);
    Gs = handleFaceDir(Gs);
    dir = cellfunUniOut(@(x)x.cells.faces(:,2), Gs);
    G.cells.faces = [cell2mat(cf), cell2mat(dir)];
    
    % geometries and layers ---------------
    geocell = {'volumes', 'centroids', 'layers'};
    for g = 1 : length(geocell)
        try
            values = cellfunUniOut(@(x)x.cells.(geocell{g}), Gs);
            G.cells.(geocell{g}) = cell2mat(values);
        catch
        end
    end

    % grdID -------------------------------
    grdID_c = arrayfun(@(x, y) x * ones(y, 1), (1:length(nc))', nc, ...
        'UniformOutput', false);
    G.cells.grdID = cell2mat(grdID_c);

    % G.faces -------------------------------------------------------------
    % num --------------------------------
    G.faces.num = sum(nf);

    % nodePos ----------------------------
    nfn = cellfunUniOut(@(x)diff(x.faces.nodePos), Gs);
    nfn = cell2mat(nfn);
    G.faces.nodePos = cumsum([1; nfn]);

    % nodes ------------------------------
    nn = cellfun(@(x)x.nodes.num, Gs);
    nn_cumsum = cumsum([0; nn]);
    nn_cumsum = mat2cell(nn_cumsum, ones(length(nn_cumsum), 1));
    fn = cellfunUniOut(@(x)x.faces.nodes(:,1), Gs);
    fn = cellfun(@(x, y) x + y, fn, nn_cumsum(1:end-1), 'UniformOutput', false);
    G.faces.nodes = cell2mat(fn);

    % neighbors --------------------------
    nc_cumsum = cumsum([0; nc]);
    neighbors = cellfunUniOut(@(x)x.faces.neighbors, Gs);
    for i = 1 : length(neighbors)
        idx = neighbors{i} ~= 0;
        neighbors{i}(idx) = neighbors{i}(idx) + nc_cumsum(i);
    end
    G.faces.neighbors = cell2mat(neighbors);

    % geometries and surfaces ---------------
    geoface = {'areas', 'normals', 'centroids', 'surfaces'};
    for g = 1 : length(geoface)
        try
            values = cellfunUniOut(@(x)x.faces.(geoface{g}), Gs);
            G.faces.(geoface{g}) = cell2mat(values);
        catch
        end
    end

    % grdID ------------------------------
    grdID_f = arrayfun(@(x, y) x * ones(y, 1), (1:length(nf))', nf, ...
        'UniformOutput', false);
    G.faces.grdID = cell2mat(grdID_f);

    % G.nodes -------------------------------------------------------------
    % num --------------------------------
    G.nodes.num = sum(nn);

    % coords -----------------------------
    coords = cellfunUniOut(@(x)x.nodes.coords, Gs);
    G.nodes.coords = cell2mat(coords);

    % grdID ------------------------------
    grdID_n = arrayfun(@(x, y) x * ones(y, 1), (1:length(nn))', nn, ...
        'UniformOutput', false);
    G.nodes.grdID = cell2mat(grdID_n);

    % G.griddim -----------------------------------------------------------

    % G.otherFields -------------------------------------------------------
    griddim = unique( cellfun(@(g)g.griddim, Gs) );
    assert(numel(griddim)==1, 'The grid dimensions are not consistent')
    G.griddim  = griddim;
    G.subGrids = Gs;

    % G.type --------------------------------------------------------------
    G.type = mfilename;
end

function Gs = handleFaceDir(Gs)
    for i = 1 : numel(Gs)
        if size(Gs{i}.cells.faces,2) < 2
            Gs{i}.cells.faces = [Gs{i}.cells.faces, ...
                nan(size(Gs{i}.cells.faces))];
        end
    end
end