function [tbls, mappings] = setupMpfaStandardTables(G, varargin)
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

    map = TensorMap();
    map.fromTbl = celltbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'cells'};
    cell_from_cellnode = map.getDispatchInd();

    map = TensorMap();
    map.fromTbl = cellnodetbl;
    map.toTbl = cellnodefacetbl;
    map.mergefds  = {'cells', 'nodes'};
    cellnode_from_cellnodeface = map.getDispatchInd();

    cell_from_cellnodeface = cell_from_cellnode(cellnode_from_cellnodeface);

    map = TensorMap();
    map.fromTbl = nodefacetbl;
    map.toTbl = cellnodefacetbl;
    map.mergefds = {'faces', 'nodes'};
    nodeface_from_cellnodeface = map.getDispatchInd();

    vectbl.vec = (1 : dim)';
    vectbl = IndexArray(vectbl);

    gen = CrossIndexArrayGenerator();
    gen.tbl1        = vectbl;
    gen.tbl2        = vectbl;
    gen.replacefds1 = {{'vec', 'vec1'}};
    gen.replacefds2 = {{'vec', 'vec2'}};
    vec12tbl = gen.eval();

    % The following tables are not virtual because they are used in scripts (there, we do not want to care about the mappings)
    cellvectbl   = crossIndexArray(celltbl, vectbl  , {}, 'optpureproduct', true);
    nodevectbl   = crossIndexArray(nodetbl, vectbl  , {}, 'optpureproduct', true);
    cellvec12tbl = crossIndexArray(celltbl, vec12tbl, {}, 'optpureproduct', true);
    
    cellnodevec12tbl   = crossIndexArray(cellnodetbl    , vec12tbl, {}, 'optpureproduct', true, 'virtual', useVirtual);
    cellnodevectbl     = crossIndexArray(cellnodetbl    , vectbl  , {}, 'optpureproduct', true, 'virtual', useVirtual);
    cellnodefacevectbl = crossIndexArray(cellnodefacetbl, vectbl  , {}, 'optpureproduct', true, 'virtual', useVirtual);
    
    cellnodeface2tbl = crossIndexArray(cellnodefacetbl, cellnodefacetbl, {'cells', 'nodes'}, 'crossextend', {{'faces', {'faces1', 'faces2'}}});

    gen = CrossIndexArrayGenerator();
    gen.tbl1        = cellnodefacetbl;
    gen.tbl2        = cellnodefacetbl;
    gen.replacefds1 = {{'faces', 'faces1'}};
    gen.replacefds2 = {{'faces', 'faces2'}};
    gen.mergefds    = {'cells', 'nodes'};
    [cellnodeface12tbl, gen] = gen.eval();

    cellnodeface1_from_cellnodeface12 = gen.ind1;
    cellnodeface2_from_cellnodeface12 = gen.ind2;
    
    cellnodeface12vectbl = crossIndexArray(cellnodeface12tbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);
     
    cellnodefacevec12tbl  = crossIndexArray(cellnodefacetbl, vec12tbl, {}, 'optpureproduct', true, 'virtual', useVirtual);

    gen = CrossIndexArrayGenerator();
    gen.tbl1        = nodefacetbl;
    gen.tbl2        = nodefacetbl;
    gen.replacefds1 = {{'faces', 'faces1'}};
    gen.replacefds2 = {{'faces', 'faces2'}};
    gen.mergefds    = { 'nodes'};
    [nodeface12tbl, gen] = gen.eval();

    nodeface1_from_nodeface12 = gen.ind1;
    nodeface2_from_nodeface12 = gen.ind2;
    
    nodefacevectbl  = crossIndexArray(nodefacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);

    tbls = struct('vec12tbl'            , vec12tbl            , ...
                  'cellvec12tbl'        , cellvec12tbl        , ...         
                  'cellnodevec12tbl'    , cellnodevec12tbl    , ...     
                  'cellnodevectbl'      , cellnodevectbl      , ...
                  'cellnodetbl'         , cellnodetbl         , ...
                  'cellnodeface12vectbl', cellnodeface12vectbl, ... 
                  'cellnodeface12tbl'   , cellnodeface12tbl   , ...    
                  'cellnodefacevec12tbl', cellnodefacevec12tbl, ... 
                  'cellnodefacevectbl'  , cellnodefacevectbl  , ...   
                  'cellnodefacetbl'     , cellnodefacetbl     , ...      
                  'celltbl'             , celltbl             , ...              
                  'facetbl'             , facetbl             , ...
                  'cellfacetbl'         , cellfacetbl         , ...
                  'vectbl'              , vectbl              , ...               
                  'nodeface12tbl'       , nodeface12tbl       , ...        
                  'nodefacevectbl'      , nodefacevectbl      , ...       
                  'nodefacetbl'         , nodefacetbl);
    
    mappings = struct('cellnode_from_cellnodeface'       , cellnode_from_cellnodeface       , ...
                      'cell_from_cellnodeface'           , cell_from_cellnodeface           , ...
                      'nodeface_from_cellnodeface'       , nodeface_from_cellnodeface       , ...
                      'cellnodeface1_from_cellnodeface12', cellnodeface1_from_cellnodeface12, ...
                      'cellnodeface2_from_cellnodeface12', cellnodeface2_from_cellnodeface12, ...
                      'nodeface1_from_nodeface12'        , nodeface1_from_nodeface12        , ...
                      'nodeface2_from_nodeface12'        , nodeface2_from_nodeface12);
    
end
