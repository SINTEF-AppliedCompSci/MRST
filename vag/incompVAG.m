function state = incompVAG(G, vagstruct, W, varargin)
% Solve incompressible flow problem (fluxes/pressures) using VAG method.
% Only Neumann bc and well input.
%
% SYNOPSIS:
%   state = incompVAG(G, vagstruct, W, varargin)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell pressures at the next time step in a
%   sequential splitting scheme for the reservoir simulation problem
%   defined by Darcy's law and a given set of external influences (wells,
%   sources, and boundary conditions).
%
%   This function uses a Vertex Approximate Gradient (VAG) method 
%
% REQUIRED PARAMETERS:
%
%   G,    - Grid and half-transmissibilities as computed by the function
%            'computeMultiPointTrans'.
%
%  vagstruct - Computed by computeVagTrans
%
% OPTIONAL PARAMETERS:
%   wells  - Well structure as defined by functions 'addWell' and
%            'assembleWellSystem'.  May be empty (i.e., W = struct([]))
%            which is interpreted as a model without any wells.
%
%
%   LinSolve     - Handle to linear system solver software to which the
%                  fully assembled system of linear equations will be
%                  passed.  Assumed to support the syntax
%
%                        x = LinSolve(A, b)
%
%                  in order to solve a system Ax=b of linear equations.
%                  Default value: LinSolve = @mldivide (backslash).
%

%
% RETURNS:

%
% SEE ALSO:
%   `computeVagTrans`, `addWell`

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opt = struct('LinSolve', @mldivide,...
                 'Verbose', mrstVerbose, ...
                 'outputFlux', false); 
    opt = merge_options(opt, varargin{:}); 

    is_well_posed = false; % changed to true if pressure is set through well
                           % or boundary conditions.
    nc = G.cells.num; 
    nn = G.nodes.num;
    
    [A, op] = setupSystem(vagstruct, G);
    
    ms = op.matrices;
    A11 = ms.A11;
    A12 = ms.A12;
    A22 = ms.A22;
    
    tbls = op.tbls;
    node2tbl    = tbls.node2tbl;
    cellnodetbl = tbls.cellnodetbl;
    cell2tbl    = tbls.cell2tbl;
    
    nw  = length(W); 
    
    welltbl.wells = (1 : nw)';
    welltbl.num = nw;
    
    wellcells = W(1).cells;
    wellind = repmat(1, numel(wellcells), 1);
    
    for iw = 2 : nw
        
        wcells = W(iw).cells;
        nwc = numel(wcells);
        wind = repmat(iw, numel(wcells), 1);
        
        wellcells = [wellcells; wcells];
        wellind   = [wellind; wind];
        
    end
    
    cellwelltbl.cells = wellcells;
    cellwelltbl.wells = wellind;
    cellwelltbl.num = numel(cellwelltbl.cells);
    
    a = A12;
    [~, cellnodewelltbl] = setupTableMapping(cellwelltbl, cellnodetbl, ...
                                                          {'cells'});
    map = setupTableMapping(cellnodetbl, cellnodewelltbl, {'cells', ...
                        'nodes'});
    a = map*a; % a belongs to cellnodewelltbl
    nodewelltbl = projTable(cellnodewelltbl, {'nodes', 'wells'});
    map = setupTableMapping(cellnodewelltbl, nodewelltbl, {'nodes', ...
                        'wells'});
    c = map*a; % c belongs to nodewelltbl;
    
    b = A22;
    [~, cellwell2tbl] = setupTableMapping(cell2tbl, cellwelltbl, {{'cells1', ...
                        'cells'}});
    cellwell2tbl = replacefield(cellwell2tbl, {'wells', 'wells1'});
    [~, cellwell2tbl] = setupTableMapping(cellwell2tbl, cellwelltbl, {{'cells2', ...
                        'cells'}});
    cellwell2tbl = replacefield(cellwell2tbl, {'wells', 'wells2'});
    map = setupTableMapping(cell2tbl, cellwell2tbl, {'cells1', 'cells2'});
    b = map*b; % b belongs to cellwell2tbl
    
    invb = 1./b;
    
    repfds = {{'cells', 'cells1'}, ...
              {'nodes', 'nodes1'}, ...
              {'wells', 'wells1'}};
    cell_1node_1well_1tbl = replacefield(cellnodewelltbl, repfds);
    repfds = {{'cells', 'cells2'}, ...
              {'nodes', 'nodes2'}, ...
              {'wells', 'wells2'}};
    cell_2node_2well_2tbl = replacefield(cellnodewelltbl, repfds);
    
    [~, prodtbl] = setupTableMapping(cellwell2tbl, cell_1node_1well_1tbl, ...
                                                   {'cells1', 'wells1'});
    [~, prodtbl] = setupTableMapping(prodtbl, cell_2node_2well_2tbl, ...
                                                   {'cells2', 'wells2'});
    map1 = setupTableMapping(cellwell2tbl, prodtbl, {'cells1', 'cells2', ...
                        'wells1', 'wells2'});
    map2 = setupTableMapping(cell_1node_1well_1tbl, prodtbl, {'cells1', 'nodes1', ...
                        'wells1'});    
    map3 = setupTableMapping(cell_2node_2well_2tbl, prodtbl, {'cells2', 'nodes2', ...
                        'wells2'});
    d = (map1*invb).*(map2*a).*(map3*a);
    e = (map1*b).*(map2*a).*(map3*a);
    
    prod_node2tbl = projTable(prodtbl, {'nodes1', 'nodes2'});
    map = setupTableMapping(prodtbl, prod_node2tbl, {'nodes1', 'nodes2'});
    d = map*d; % d belongs to prod_node2tbl
    prod_well2tbl = projTable(prodtbl, {'wells1', 'wells2'});
    map = setupTableMapping(prodtbl, prod_well2tbl, {'wells1', 'wells2'});
    e = map*e; % e belongs to prod_well2tbl
    

    tbl = prod_node2tbl; % alias
    Ann = sparse(tbl.nodes1, tbl.nodes2, d, nn, nn);
    tbl = nodewelltbl; % alias    
    Anw = sparse(tbl.nodes, tbl.wells, c, nn, nw);
    Awn = Anw';
    tbl = prod_well2tbl; % alias
    Aww = sparse(tbl.wells1, tbl.wells2, e, nw, nw);
    
    Awq = speye(nw);
    Anq = sparse(nn, nw);
    Aext = [[A + Ann, Anw,  Anq];
            [    Awn, Aww, -Awq];
           ];
    

    state = [];
    
end


