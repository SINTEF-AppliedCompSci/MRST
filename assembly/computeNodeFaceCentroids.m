function [cellnodefacecents, nodefacecents] = computeNodeFaceCentroids(G, eta, tbls, varargin)
    
    opt = struct('bcetazero', false);
    opt = merge_options(opt, varargin{:});
    bcetazero = opt.bcetazero;

    % cellnodefacecents centroid of node-face points, relative to cell
    % centroids. It belongs to cellnodefacecoltbl
    %
    % nodefacecents centroid of node-face points, belongs to nodefacecoltbl

    cellnodefacetbl    = tbls.cellnodefacetbl;
    cellnodefacecoltbl = tbls.cellnodefacecoltbl;
    nodefacecoltbl     = tbls.nodefacecoltbl;
    
    cno = cellnodefacetbl.get('cells');
    fno = cellnodefacetbl.get('faces');
    nno = cellnodefacetbl.get('nodes');
    
    % Absolute position of node-face points (in cellnodefacecoltbl)
    ccents = G.cells.centroids(cno, :);
    fcents = G.faces.centroids(fno, :);
    ncents = G.nodes.coords(nno, :);

    if bcetazero
        % zero eta at boundary faces
        eta = eta*ones(cellnodefacetbl.num, 1);
        bcfnoind = any(G.faces.neighbors(fno, :) == 0, 2);
        eta(bcfnoind) = 0;
    end
    
    abscellnodefacecents = bsxfun(@times, eta, ncents) + bsxfun(@times, (1 - ...
                                                      eta), fcents);
    
    % Relative position of node-face points (in cellnodefacecoltbl)
    cellnodefacecents = abscellnodefacecents - ccents;
    
    cellnodefacecents = reshape(cellnodefacecents', [], 1);
    abscellnodefacecents = reshape(abscellnodefacecents', [], 1);

    map = TensorMap();
    map.fromTbl = cellnodefacecoltbl;
    map.toTbl = nodefacecoltbl;
    map.mergefds = {'nodes', 'faces', 'coldim'};
    map = map.setup();
    
    nodefacecents = map.eval(abscellnodefacecents);
    coef = map.eval(ones(cellnodefacecoltbl.num, 1));
    
    nodefacecents = 1./coef.*nodefacecents;

end
