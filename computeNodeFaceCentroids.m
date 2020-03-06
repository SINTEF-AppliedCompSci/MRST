function [cellnodefacecents, nodefacecents] = computeNodeFaceCentroids(G, tbls, eta)

% cellnodefacecents centroid of node-face points, relative to cell
% centroids. It belongs to cellnodefacecoltbl
%
% nodefacecents centroid of node-face points, belongs to nodefacecoltbl
    
    cellnodefacetbl    = tbls.cellnodefacetbl;
    cellnodefacecoltbl = tbls.cellnodefacecoltbl;
    nodefacecoltbl     = tbls.nodefacecoltbl;
    
    fno = cellnodefacetbl.get('faces');
    cno = cellnodefacetbl.get('cells');
    nno = cellnodefacetbl.get('nodes');
    
    % Absolute position of node-face points (in cellnodefacecoltbl)
    abscellnodefacecents = G.faces.centroids(fno, :) + eta* (G.nodes.coords(nno, ...
                                                      :) - G.faces.centroids(fno, ...
                                                      :));
    % Relative position of node-face points (in cellnodefacecoltbl)
    cellnodefacecents = abscellnodefacecents - G.cells.centroids(cno, :);
    
    
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
