function [tbls, mappings] = setupStandardTables(G, varargin)
    
    opt = struct('useVirtual', true);
    opt = merge_options(opt, varargin{:});
    useVirtual = opt.useVirtual;
    
    nc  = G.cells.num;
    nf  = G.faces.num;
    nn  = G.nodes.num;
    dim = G.griddim;

    coltbl.coldim = (1 : dim)';
    coltbl = IndexArray(coltbl);
    rowtbl = coltbl;
    rowtbl = replacefield(rowtbl, {'coldim', 'rowdim'});

    celltbl.cells = (1 : nc)';
    celltbl = IndexArray(celltbl);
    
    nodetbl.nodes = (1 : nn)';
    nodetbl = IndexArray(nodetbl);
    
    cellcoltbl = crossIndexArray(celltbl, coltbl, {}); % ordering is cell - col
    nodecoltbl = crossIndexArray(nodetbl, coltbl, {}); % ordering is cell - col

    cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos)); 
    cellfacetbl.faces = G.cells.faces(:, 1);
    cellfacetbl = IndexArray(cellfacetbl);
    
    nodefacetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos)); 
    nodefacetbl.nodes = G.faces.nodes;
    nodefacetbl = IndexArray(nodefacetbl); 
    
    % We setup the face-node table and it is ordered along ascending node numbers so
    % that we will have a block structure for the nodal scalar product.
    nodefacetbl = sortIndexArray(nodefacetbl, {'nodes', 'faces'});
    
    % not virtual because used in setupBCpercase (could be optimized)
    nodefacecoltbl = crossIndexArray(nodefacetbl, coltbl, {});

    % We setup the cell-face-node table, cellnodefacetbl. Each entry determine a
    % unique facet in a corner
    % We order cellnodeface in cell-node-face order. This is node to optimize
    % for-end loop below.
    cellnodefacetbl = crossIndexArray(cellfacetbl, nodefacetbl, {'faces'});
    cellnodefacetbl = sortIndexArray(cellnodefacetbl, {'cells', 'nodes', 'faces'});

    % We setup the cell-node table, cellnodetbl. Each entry determine a unique
    % corner
    cellnodetbl = projIndexArray(cellnodefacetbl, {'nodes', 'cells'});
    cellnodetbl = sortIndexArray(cellnodetbl, {'cells', 'nodes'});

    map = TensorMap();
    map.fromTbl = celltbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'cells'};
    cell_from_cellnode = getDispatchInd(map);
    
    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'nodes'};
    node_from_cellnode = getDispatchInd(map);
    
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = cellnodefacetbl;
    map.mergefds  = {'cells', 'nodes'};
    cellnode_from_cellnodeface = getDispatchInd(map);
    
    map = TensorMap();
    map.fromTbl = nodefacetbl;
    map.toTbl = cellnodefacetbl;
    map.mergefds = {'faces', 'nodes'};
    nodeface_from_cellnodeface = getDispatchInd(map);

    cellnodecoltbl    = crossIndexArray(cellnodetbl, coltbl, {}, 'virtual', useVirtual);

    % not virtual because used in setupBCpercase (could be optimized)
    cellnodefacecoltbl = crossIndexArray(cellnodefacetbl, coltbl, {});
    
    colrowtbl = crossIndexArray(coltbl, rowtbl, {});
    cellnodecolrowtbl = crossIndexArray(cellnodetbl, colrowtbl, {}, 'virtual', ...
                                        useVirtual);
    cellnodefacecolrowtbl = crossIndexArray(cellnodefacetbl, colrowtbl, {}, ...
                                            'virtual', useVirtual);

    nodecolrowtbl = crossIndexArray(nodetbl, colrowtbl, {}, 'virtual', useVirtual);
    
    fds = {{'rowdim', {'rowdim1', 'rowdim2'}}, ...
           {'coldim', {'coldim1', 'coldim2'}}};
    col2row2tbl = crossIndexArray(colrowtbl, colrowtbl, {}, 'crossextend', fds);
    
    % not virtual because used in setupStiffnessTensor (could be optimized).
    cellcol2row2tbl = crossIndexArray(celltbl, col2row2tbl, {}, 'optpureproduct', ...
                                      true);
    cellnodecol2row2tbl = crossIndexArray(cellnodetbl, col2row2tbl, {}, ...
                                          'virtual', useVirtual);
    
    tbls = struct('coltbl'               , coltbl               , ...
                  'celltbl'              , celltbl              , ...
                  'nodetbl'              , nodetbl              , ...
                  'cellfacetbl'          , cellfacetbl          , ...
                  'cellnodetbl'          , cellnodetbl          , ...
                  'nodefacetbl'          , nodefacetbl          , ...
                  'cellcoltbl'           , cellcoltbl           , ... 
                  'nodecoltbl'           , nodecoltbl           , ... 
                  'nodefacecoltbl'       , nodefacecoltbl       , ... 
                  'cellnodefacetbl'      , cellnodefacetbl      , ... 
                  'cellnodecoltbl'       , cellnodecoltbl       , ...    
                  'cellnodecolrowtbl'    , cellnodecolrowtbl    , ... 
                  'cellnodefacecoltbl'   , cellnodefacecoltbl   , ... 
                  'cellnodefacecolrowtbl', cellnodefacecolrowtbl, ... 
                  'colrowtbl'            , colrowtbl            , ... 
                  'nodecolrowtbl'        , nodecolrowtbl        , ... 
                  'col2row2tbl'          , col2row2tbl          , ... 
                  'cellcol2row2tbl'      , cellcol2row2tbl      , ...
                  'cellnodecol2row2tbl'  , cellnodecol2row2tbl);

    mappings = struct('cell_from_cellnode'        , cell_from_cellnode        , ...
                      'node_from_cellnode'        , node_from_cellnode        , ...
                      'cellnode_from_cellnodeface', cellnode_from_cellnodeface, ...
                      'nodeface_from_cellnodeface', nodeface_from_cellnodeface);
    
end
