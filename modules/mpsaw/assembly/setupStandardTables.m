function [tbls, mappings] = setupStandardTables(G, varargin)
%Undocumented Utility Function

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


    opt = struct('useVirtual', true, ...
                 'inittbls', []);
    opt = merge_options(opt, varargin{:});
    useVirtual = opt.useVirtual;

    nc  = G.cells.num;
    nf  = G.faces.num;
    nn  = G.nodes.num;
    dim = G.griddim;

    if ~isempty(opt.inittbls)
        itbls = opt.inittbls;
        celltbl         = itbls.celltbl;
        facetbl         = itbls.facetbl;
        nodetbl         = itbls.nodetbl;
        cellnodetbl     = itbls.cellnodetbl;
        nodefacetbl     = itbls.nodefacetbl;
        cellfacetbl     = itbls.cellfacetbl;
        cellnodefacetbl = itbls.cellnodefacetbl;

    else
        celltbl.cells = (1 : nc)';
        celltbl = IndexArray(celltbl);
    
        facetbl.faces = (1 : nf)';
        facetbl = IndexArray(facetbl);

        nodetbl.nodes = (1 : nn)';
        nodetbl = IndexArray(nodetbl);
    
        cellfacetbl.cells = rldecode((1 : nc)', diff(G.cells.facePos)); 
        cellfacetbl.faces = G.cells.faces(:, 1);
        cellfacetbl = IndexArray(cellfacetbl);
        
        nodefacetbl.faces = rldecode((1 : nf)', diff(G.faces.nodePos)); 
        nodefacetbl.nodes = G.faces.nodes;
        nodefacetbl = IndexArray(nodefacetbl); 
    
        % We setup the face-node table and it is ordered along ascending node numbers so
        % that we will have a block structure for the nodal scalar product.
        nodefacetbl = sortIndexArray(nodefacetbl, {'nodes', 'faces'});
    
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
    
    end
    
    coltbl.coldim = (1 : dim)';
    coltbl = IndexArray(coltbl);
    rowtbl = coltbl;
    rowtbl = replacefield(rowtbl, {'coldim', 'rowdim'});

    cellcoltbl = crossIndexArray(celltbl, coltbl, {}); % ordering is cell - col
    nodecoltbl = crossIndexArray(nodetbl, coltbl, {}); % ordering is node - col
    % not virtual because used in setupBCpercase (could be optimized)
    nodefacecoltbl = crossIndexArray(nodefacetbl, coltbl, {});

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

    cell_from_cellnodeface = cell_from_cellnode(cellnode_from_cellnodeface);

    % for mpfa
    cellnodeface2tbl = crossIndexArray(cellnodefacetbl, cellnodefacetbl, {'cells', 'nodes'}, 'crossextend', {{'faces', {'faces1', 'faces2'}}});
    cellnodeface2coltbl = crossIndexArray(cellnodeface2tbl, coltbl, {}, 'virtual', useVirtual);
    
    map = TensorMap();
    map.fromTbl = cellnodefacetbl;
    map.toTbl = cellnodeface2tbl;
    map.replaceFromTblfds = {{'faces', 'faces1'}};
    map.mergefds = {'cells', 'nodes', 'faces1'};
    
    cellnodeface_1_from_cellnodeface2 = map.getDispatchInd();
    
    map = TensorMap();
    map.fromTbl = cellnodefacetbl;
    map.toTbl = cellnodeface2tbl;
    map.replaceFromTblfds = {{'faces', 'faces2'}};
    map.mergefds = {'cells', 'nodes', 'faces2'};
    
    cellnodeface_2_from_cellnodeface2 = map.getDispatchInd();
    
    nodeface2tbl = crossIndexArray(nodefacetbl, nodefacetbl, {'nodes'}, 'crossextend', {{'faces', {'faces1', 'faces2'}}});
    
    map = TensorMap();
    map.fromTbl = nodefacetbl;
    map.toTbl = nodeface2tbl;
    map.replaceFromTblfds = {{'faces', 'faces1'}};
    map.mergefds = {'nodes', 'faces1'};
    
    nodeface_1_from_nodeface2 = map.getDispatchInd();
    
    map = TensorMap();
    map.fromTbl = nodefacetbl;
    map.toTbl = nodeface2tbl;
    map.replaceFromTblfds = {{'faces', 'faces2'}};
    map.mergefds = {'nodes', 'faces2'};
    
    nodeface_2_from_nodeface2 = map.getDispatchInd();
    
    
    cellnodecoltbl = crossIndexArray(cellnodetbl, coltbl, {}, 'virtual', useVirtual);

    % not virtual because used in setupBCpercase (could be optimized)
    cellnodefacecoltbl = crossIndexArray(cellnodefacetbl, coltbl, {});
    
    colrowtbl = crossIndexArray(coltbl, rowtbl, {});
    
    cellcolrowtbl = crossIndexArray(celltbl, colrowtbl, {}, 'virtual', ...
                                    useVirtual);
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
                  'rowtbl'               , rowtbl               , ... 
                  'celltbl'              , celltbl              , ...
                  'facetbl'              , facetbl              , ...
                  'nodetbl'              , nodetbl              , ...
                  'cellfacetbl'          , cellfacetbl          , ...
                  'cellnodetbl'          , cellnodetbl          , ...
                  'nodefacetbl'          , nodefacetbl          , ...
                  'cellcoltbl'           , cellcoltbl           , ... 
                  'cellcolrowtbl'        , cellcolrowtbl        , ...
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
                  'cellnodecol2row2tbl'  , cellnodecol2row2tbl  , ....
                  'cellnodeface2tbl'     , cellnodeface2tbl     , ...
                  'cellnodeface2coltbl'  , cellnodeface2coltbl     , ...
                  'nodeface2tbl'         , nodeface2tbl);
    
    mappings = struct('cell_from_cellnode'               , cell_from_cellnode               , ...
                      'node_from_cellnode'               , node_from_cellnode               , ...
                      'cell_from_cellnodeface'           , cell_from_cellnodeface           , ...
                      'cellnode_from_cellnodeface'       , cellnode_from_cellnodeface       , ...
                      'nodeface_from_cellnodeface'       , nodeface_from_cellnodeface       , ...
                      'cellnodeface_1_from_cellnodeface2', cellnodeface_1_from_cellnodeface2, ... 
                      'cellnodeface_2_from_cellnodeface2', cellnodeface_2_from_cellnodeface2, ... 
                      'nodeface_1_from_nodeface2'        , nodeface_1_from_nodeface2        , ...         
                      'nodeface_2_from_nodeface2'        , nodeface_2_from_nodeface2);
    
    
    tbls.useVirtual = useVirtual;
    
end
