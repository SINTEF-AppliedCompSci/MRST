function [tbls, mappings] = setupStandardBlockTables(G, nodetbl, globtbls, varargin)
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
                 'tblcase', 'empty');
    
    opt = merge_options(opt, varargin{:});
    useVirtual = opt.useVirtual;

    globcelltbl         = globtbls.celltbl;
    globfacetbl         = globtbls.facetbl;
    globnodetbl         = globtbls.nodetbl;
    globcellnodetbl     = globtbls.cellnodetbl;
    globnodefacetbl     = globtbls.nodefacetbl;
    globcellnodefacetbl = globtbls.cellnodefacetbl;
    
    cellnodetbl = crossIndexArray(nodetbl, globcellnodetbl, {'nodes'});
    cellnodetbl = sortIndexArray(cellnodetbl, {'cells', 'nodes'});
    
    celltbl = projIndexArray(cellnodetbl, {'cells'});
    
    map = TensorMap();
    map.fromTbl  = globcelltbl;
    map.toTbl    = celltbl;    
    map.mergefds = {'cells'};
    globcell_from_cell = map.getDispatchInd();
    
    
    nodefacetbl = crossIndexArray(nodetbl, globnodefacetbl, {'nodes'});
    nodefacetbl = sortIndexArray(nodefacetbl, {'nodes', 'faces'});    
    
    facetbl = projIndexArray(nodefacetbl, {'faces'});
    
    cellnodefacetbl = crossIndexArray(nodetbl, globcellnodefacetbl, {'nodes'});
    cellnodefacetbl = sortIndexArray(cellnodefacetbl, {'cells', 'nodes', ...
                        'faces'});
    
    cellfacetbl = projIndexArray(cellnodefacetbl, {'cells', 'faces'});
    
    inittbls = struct('celltbl'        , celltbl    , ...
                      'nodetbl'        , nodetbl    , ...
                      'facetbl'        , facetbl    , ...
                      'cellnodetbl'    , cellnodetbl, ...
                      'nodefacetbl'    , nodefacetbl, ...
                      'cellfacetbl'    , cellfacetbl, ...
                      'cellnodefacetbl', cellnodefacetbl);

    switch opt.tblcase
      case 'mpsa'
        [tbls, mappings] = setupMpsaStandardTables(G, 'inittbls', inittbls, 'useVirtual', useVirtual);
      case 'mpfa'
        [tbls, mappings] = setupMpfaStandardTables(G, 'inittbls', inittbls, 'useVirtual', useVirtual);
      case 'both'
        [tbls, mappings] = setupMpxaStandardTables(G, 'inittbls', inittbls, 'useVirtual', useVirtual);
      case 'empty'
        error('to avoid unnecessary index array setup, provide the tblcase option [mpsa/mpfa/both]')
      otherwise
        error('tblcase not recognized')
    end

    mappings.globcell_from_cell = globcell_from_cell;

    % add global tables and mappings in output
    if useVirtual
        
        map = TensorMap();
        map.fromTbl  = globnodefacetbl;
        map.toTbl    = nodefacetbl;
        map.mergefds = {'nodes', 'faces'};

        mappings.globnodeface_from_nodeface = map.getDispatchInd();

        map = TensorMap();
        map.fromTbl  = globfacetbl;
        map.toTbl    = facetbl;    
        map.mergefds = {'faces'};
        mappings.globface_from_face = map.getDispatchInd();
    
        
    end
    
    tbls.globcelltbl = globcelltbl;

end
