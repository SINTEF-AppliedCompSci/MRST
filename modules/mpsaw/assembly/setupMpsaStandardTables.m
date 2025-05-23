function [tbls, mappings] = setupMpsaStandardTables(G, varargin)
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


    opt = struct('useVirtual', false, ...
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
    
    vectbl.vec = (1 : dim)';
    vectbl = IndexArray(vectbl);

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = vectbl;
    gen.tbl2 = vectbl;
    gen.replacefds1 = {{'vec', 'vec1'}};
    gen.replacefds2 = {{'vec', 'vec2'}};

    vec12tbl = gen.eval();

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = vectbl;
    gen.tbl2 = vec12tbl;
    gen.replacefds1 = {{'vec', 'vec3'}};

    vec123tbl = gen.eval();
    vec123tbl = sortIndexArray(vec123tbl, {'vec1', 'vec2', 'vec3'});
    
    map = TensorMap();
    map.fromTbl = celltbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'cells'};
    cell_from_cellnode = map.getDispatchInd();

    map = TensorMap();
    map.fromTbl = facetbl;
    map.toTbl = cellfacetbl;
    map.mergefds = {'faces'};
    face_from_cellface = map.getDispatchInd();

    map = TensorMap();
    map.fromTbl = celltbl;
    map.toTbl = cellfacetbl;
    map.mergefds = {'cells'};
    cell_from_cellface = map.getDispatchInd();
    
    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'nodes'};
    node_from_cellnode = map.getDispatchInd();

    map = TensorMap();
    map.fromTbl = facetbl;
    map.toTbl = nodefacetbl;
    map.mergefds = {'faces'};
    face_from_nodeface = map.getDispatchInd();

    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = nodefacetbl;
    map.mergefds = {'nodes'};
    node_from_nodeface = map.getDispatchInd();
    
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
    
    cellnodevectbl     = crossIndexArray(cellnodetbl    , vectbl  , {}, 'optpureproduct', true, 'virtual', useVirtual);
    cellvectbl         = crossIndexArray(celltbl        , vectbl  , {}, 'optpureproduct', true, 'virtual', useVirtual); 
    nodevectbl         = crossIndexArray(nodetbl        , vectbl  , {}, 'optpureproduct', true, 'virtual', useVirtual); 
    cellnodefacevectbl = crossIndexArray(cellnodefacetbl, vectbl  , {}, 'optpureproduct', true, 'virtual', useVirtual);
    nodefacevectbl     = crossIndexArray(nodefacetbl    , vectbl  , {}, 'optpureproduct', true, 'virtual', useVirtual);
    facevectbl         = crossIndexArray(facetbl        , vectbl  , {}, 'optpureproduct', true, 'virtual', useVirtual);
    cellnodevec12tbl   = crossIndexArray(cellnodetbl    , vec12tbl, {}, 'optpureproduct', true, 'virtual', useVirtual);
    cellfacevectbl     = crossIndexArray(cellfacetbl    , vectbl  , {}, 'optpureproduct', true, 'virtual', useVirtual);
    cellvec12tbl       = crossIndexArray(celltbl        , vec12tbl, {}, 'optpureproduct', true, 'virtual', useVirtual);
    
    gen = CrossIndexArrayGenerator();
    gen.tbl1 = vec12tbl;
    gen.tbl2 = vec12tbl;
    gen.replacefds1 = {{'vec1', 'vec11'},{'vec2', 'vec12'}};
    gen.replacefds2 = {{'vec1', 'vec21'},{'vec2', 'vec22'}};

    vec1212tbl = gen.eval();

    % not virtual because used in setupStiffnessTensor.
    cellvec1212tbl = crossIndexArray(celltbl, vec1212tbl, {}, 'optpureproduct', true, 'virtual', useVirtual);

    gen = CrossIndexArrayGenerator();
    gen.tbl1        = cellnodetbl;
    gen.tbl2        = cellnodetbl;
    gen.replacefds1 = {{'cells', 'cells1'}};
    gen.replacefds2 = {{'cells', 'cells2'}};
    gen.mergefds    = {'nodes'};

    [cell12nodetbl, gen] = gen.eval();

    cell1node_from_cell12node = gen.ind1;
    cell2node_from_cell12node = gen.ind2;
    
    gen = CrossIndexArrayGenerator();
    gen.tbl1        = vec12tbl;
    gen.tbl2        = vectbl;
    gen.replacefds1 = {{'vec1', 'vec11'}, {'vec2', 'vec12'}};
    gen.replacefds2 = {{'vec', 'vec2'}};
    gen.opts        = {'optpureproduct', true};
    vec122tbl = gen.eval();
    
    cellnodefacevec122tbl = crossIndexArray(cellnodefacetbl, vec122tbl, {}, ...
                                            'optpureproduct', true, ...
                                            'virtual', useVirtual);

    gen = CrossIndexArrayGenerator();
    gen.tbl1        = cellnodetbl;
    gen.tbl2        = cellnodefacetbl;
    gen.replacefds1 = {{'cells', 'cells1'}};
    gen.replacefds2 = {{'cells', 'cells2'}};
    gen.mergefds    = {'nodes'};
    
    [cell12nodefacetbl, gen] = gen.eval();
    cell1node_from_cell12nodeface     = gen.ind1;
    cell2nodeface_from_cell12nodeface = gen.ind2;
    
    cell12nodefacevec122tbl = crossIndexArray(cell12nodefacetbl, vec122tbl, {}, ...
                                              'optpureproduct', true, ...
                                              'virtual', useVirtual);

    map = TensorMap();
    map.fromTbl  = cell12nodetbl;
    map.toTbl    = cell12nodefacetbl;
    map.mergefds = {'cells1', 'cells2', 'nodes'};

    cell12node_from_cell12nodeface = map.getDispatchInd();
    
    gen = CrossIndexArrayGenerator();
    gen.tbl1        = cellnodefacetbl; 
    gen.tbl2        = cell12nodefacetbl;
    gen.replacefds1 = {{'faces', 'faces1'}, {'cells', 'cells1'}};
    gen.replacefds2 = {{'faces', 'faces2'}};
    gen.mergefds    = {'cells1', 'nodes'};
    
    [cell12nodeface12tbl, gen] = gen.eval();
    cell1nodeface1_from_cell12nodeface12  = gen.ind1;
    cell12nodeface2_from_cell12nodeface12 = gen.ind2;
    
    cell12nodeface12vec12tbl = crossIndexArray(cell12nodeface12tbl, vec12tbl, {}, ...
                                               'optpureproduct', true, ...
                                               'virtual', useVirtual);
    
    tbls = struct('vectbl'                  , vectbl                  , ...
                  'celltbl'                 , celltbl                 , ...
                  'facetbl'                 , facetbl                 , ...
                  'nodetbl'                 , nodetbl                 , ...
                  'cellfacetbl'             , cellfacetbl             , ...
                  'cellnodetbl'             , cellnodetbl             , ...
                  'nodefacetbl'             , nodefacetbl             , ...
                  'cellvectbl'              , cellvectbl              , ... 
                  'cellvec12tbl'            , cellvec12tbl            , ...
                  'nodevectbl'              , nodevectbl              , ... 
                  'nodefacevectbl'          , nodefacevectbl          , ... 
                  'cellnodefacetbl'         , cellnodefacetbl         , ... 
                  'cellnodevectbl'          , cellnodevectbl          , ...    
                  'cellnodevec12tbl'        , cellnodevec12tbl        , ... 
                  'cellnodefacevectbl'      , cellnodefacevectbl      , ... 
                  'vec12tbl'                , vec12tbl                , ... 
                  'vec123tbl'               , vec123tbl               , ...
                  'vec122tbl'               , vec122tbl               , ...
                  'vec1212tbl'              , vec1212tbl              , ...
                  'facevectbl'              , facevectbl              , ...
                  'cellfacevectbl'          , cellfacevectbl          , ...
                  'cell12nodetbl'           , cell12nodetbl           , ...
                  'cellnodefacevec122tbl'   , cellnodefacevec122tbl   , ...
                  'cell12nodefacetbl'       , cell12nodefacetbl       , ...
                  'cell12nodeface12tbl'     , cell12nodeface12tbl     , ...
                  'cell12nodefacevec122tbl' , cell12nodefacevec122tbl , ...
                  'cell12nodeface12vec12tbl', cell12nodeface12vec12tbl, ...
                  'cellvec1212tbl'          , cellvec1212tbl);
    
    mappings = struct('cell_from_cellnode'                   , cell_from_cellnode                   , ...
                      'face_from_cellface'                   , face_from_cellface                   , ...
                      'cell_from_cellface'                   , cell_from_cellface                   , ...
                      'face_from_nodeface'                   , face_from_nodeface                   , ...
                      'node_from_nodeface'                   , node_from_nodeface                   , ...
                      'node_from_cellnode'                   , node_from_cellnode                   , ...
                      'cell_from_cellnodeface'               , cell_from_cellnodeface               , ...
                      'cellnode_from_cellnodeface'           , cellnode_from_cellnodeface           , ...
                      'nodeface_from_cellnodeface'           , nodeface_from_cellnodeface           , ...
                      'cell1node_from_cell12nodeface'        , cell1node_from_cell12nodeface        , ...
                      'cell2nodeface_from_cell12nodeface'    , cell2nodeface_from_cell12nodeface    , ...
                      'cell12node_from_cell12nodeface'       , cell12node_from_cell12nodeface       , ...
                      'cell1nodeface1_from_cell12nodeface12' , cell1nodeface1_from_cell12nodeface12 , ... 
                      'cell12nodeface2_from_cell12nodeface12', cell12nodeface2_from_cell12nodeface12, ...
                      'cell1node_from_cell12node'            , cell1node_from_cell12node            , ...
                      'cell2node_from_cell12node'            , cell2node_from_cell12node);
    
    tbls.useVirtual = useVirtual;
    
end
