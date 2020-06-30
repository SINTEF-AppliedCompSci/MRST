function [tbls, mappings] = setupStandardBlockTables(G, nodetbl, globtbls, varargin)
    
    opt = struct('useVirtual', true);
    opt = merge_options(opt, varargin{:});
    useVirtual = opt.useVirtual;

    globcelltbl         = globtbls.celltbl;
    globnodetbl         = globtbls.nodetbl;
    globcellnodetbl     = globtbls.cellnodetbl;
    globnodefacetbl     = globtbls.nodefacetbl;
    globcellnodefacetbl = globtbls.cellnodefacetbl;
    
    cellnodetbl = crossIndexArray(nodetbl, globcellnodetbl, {'nodes'});
    cellnodetbl = sortIndexArray(cellnodetbl, {'cells', 'nodes'});
    
    celltbl = projIndexArray(cellnodetbl, {'cells'});
    
    map = TensorMap();
    map.fromTbl = globcelltbl;
    map.toTbl = celltbl;    
    map.mergefds = {'cells'};
    globcell_from_cell = map.getDispatchInd();
    
    
    nodefacetbl = crossIndexArray(nodetbl, globnodefacetbl, {'nodes'});
    nodefacetbl = sortIndexArray(nodefacetbl, {'nodes', 'faces'});    
    
    facetbl = projIndexArray(nodefacetbl, {'faces'});
    
    cellnodefacetbl = crossIndexArray(nodetbl, globcellnodefacetbl, {'nodes'});
    cellnodefacetbl = sortIndexArray(cellnodefacetbl, {'cells', 'nodes', ...
                        'faces'});
    
    cellfacetbl = projIndexArray(cellnodefacetbl, {'cells', 'faces'});
    
    inittbls = struct('celltbl', celltbl, ...
                      'nodetbl', nodetbl, ...
                      'facetbl', facetbl, ...
                      'cellnodetbl', cellnodetbl, ...
                      'nodefacetbl', nodefacetbl, ...
                      'cellfacetbl', cellfacetbl, ...
                      'cellnodefacetbl', cellnodefacetbl);
        
    [tbls, mappings] = setupStandardTables(G, 'inittbls', inittbls, 'useVirtual', useVirtual);
    
    mappings.globcell_from_cell = globcell_from_cell;
        
end
