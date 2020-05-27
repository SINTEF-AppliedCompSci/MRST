function [tbls, mappings] = setupStandardBlockTables(G, nodetbl, globtbls, varargin)
    
    opt = struct('useVirtual', true);
    opt = merge_options(opt, varargin{:});
    useVirtual = opt.useVirtual;

    globcelltbl         = globtbls.celltbl;
    globnodetbl         = globtbls.nodetbl;
    globcellnodetbl     = globtbls.cellnodetbl;
    globnodefacetbl     = globtbls.nodefacetbl;
    globcellnodefacetbl = globtbls.cellnodefacetbl;
    coltbl              = globtbls.coltbl;
    colrowtbl           = globtbls.colrowtbl;
    col2row2tbl         = globtbls.col2row2tbl;
    
    cellnodetbl = crossIndexArray(nodetbl, globcellnodetbl, {'nodes'});
    cellnodetbl = sortIndexArray(cellnodetbl, {'cells', 'nodes'});
    
    celltbl = projIndexArray(cellnodetbl, {'cells'});

    map = TensorMap();
    map.fromTbl = globcelltbl;
    map.toTbl = celltbl;    
    map.mergefds = {'cells'};
    globcell_from_cell = map.getDispatchInd();
    
    map = TensorMap();
    map.fromTbl = celltbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'cells'};
    cell_from_cellnode = map.getDispatchInd();
    
    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'nodes'};
    node_from_cellnode = map.getDispatchInd(); 
   
    
    nodefacetbl = crossIndexArray(nodetbl, globnodefacetbl, {'nodes'});
    nodefacetbl = sortIndexArray(nodefacetbl, {'nodes', 'faces'});    
    
    cellnodefacetbl = crossIndexArray(nodetbl, globcellnodefacetbl, {'nodes'});
    cellnodefacetbl = sortIndexArray(cellnodefacetbl, {'cells', 'nodes', ...
                        'faces'});
    
    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = cellnodefacetbl;
    map.mergefds  = {'cells', 'nodes'};
    cellnode_from_cellnodeface = map.getDispatchInd();
    
    map = TensorMap();
    map.fromTbl = nodefacetbl;
    map.toTbl = cellnodefacetbl;
    map.mergefds = {'faces', 'nodes'};
    nodeface_from_cellnodeface = map.getDispatchInd();

    cellcoltbl = crossIndexArray(celltbl, coltbl, {}); % ordering is cell - col
    nodecoltbl = crossIndexArray(nodetbl, coltbl, {}); % ordering is cell - col

    
    % not virtual because used in setupBCpercase (could be optimized)
    nodefacecoltbl = crossIndexArray(nodefacetbl, coltbl, {});

    cellnodecoltbl    = crossIndexArray(cellnodetbl, coltbl, {}, 'virtual', useVirtual);

    % not virtual because used in setupBCpercase (could be optimized)
    cellnodefacecoltbl = crossIndexArray(cellnodefacetbl, coltbl, {});
    
    cellnodecolrowtbl = crossIndexArray(cellnodetbl, colrowtbl, {}, 'virtual', ...
                                        useVirtual);
    cellnodefacecolrowtbl = crossIndexArray(cellnodefacetbl, colrowtbl, {}, ...
                                            'virtual', useVirtual);

    nodecolrowtbl = crossIndexArray(nodetbl, colrowtbl, {}, 'virtual', useVirtual);
    cellcolrowtbl = crossIndexArray(celltbl, colrowtbl, {}, 'virtual', useVirtual);
    
    % not virtual because used in setupStiffnessTensor (could be optimized).
    cellcol2row2tbl = crossIndexArray(celltbl, col2row2tbl, {}, 'optpureproduct', ...
                                      true);
    cellnodecol2row2tbl = crossIndexArray(cellnodetbl, col2row2tbl, {}, ...
                                          'virtual', useVirtual);
    
    tbls = struct('coltbl'               , coltbl               , ...
                  'celltbl'              , celltbl              , ...
                  'nodetbl'              , nodetbl              , ...
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
                  'cellcolrowtbl'        , cellcolrowtbl        , ...
                  'nodecolrowtbl'        , nodecolrowtbl        , ... 
                  'col2row2tbl'          , col2row2tbl          , ... 
                  'cellcol2row2tbl'      , cellcol2row2tbl      , ...
                  'cellnodecol2row2tbl'  , cellnodecol2row2tbl);

    mappings = struct('cell_from_cellnode'        , cell_from_cellnode        , ...
                      'globcell_from_cell'        , globcell_from_cell        , ...
                      'node_from_cellnode'        , node_from_cellnode        , ...
                      'cellnode_from_cellnodeface', cellnode_from_cellnodeface, ...
                      'nodeface_from_cellnodeface', nodeface_from_cellnodeface);
        
end
