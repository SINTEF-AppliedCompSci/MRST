function pn = incompVAG(G, vagstruct, W, varargin)
%
% SYNOPSIS:
%   state = incompVAG(G, vagstruct, W, varargin)
%
% DESCRIPTION:
%   Solve incompressible flow problem (fluxes/pressures) using VAG method.
%   Only Neumann bc and well input.
%
% REQUIRED PARAMETERS:
%
%  G         - Grid
%  vagstruct - Structure as computed by computeVagTrans
%
% OPTIONAL PARAMETERS:
%   wells  - Well structure as defined by functions 'addWell' and
%            'assembleWellSystem'.  May be empty (i.e., W = struct([]))
%            which is interpreted as a model without any wells.
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
% RETURNS: pressure at nodes

%
% SEE ALSO:
%   `computeVagTrans`, `addWell`

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
    wellind   = repmat(1, numel(wellcells), 1);
    welltypes = getWellControlType(W(1));
    wellvals  = W(1).val;
    
    for iw = 2 : nw
        
        wcells = W(iw).cells;
        nwc = numel(wcells);
        wind = repmat(iw, numel(wcells), 1);
        wtype = getWellControlType(W(iw));
        wval = W(iw).val;
        
        wellcells = [wellcells; wcells];
        wellind   = [wellind; wind];
        welltypes = [welltypes; wtype];
        wellvals  = [wellvals; wval];
        
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
    [~, cell2well2tbl] = setupTableMapping(cell2tbl, cellwelltbl, {{'cells1', ...
                        'cells'}});
    cell2well2tbl = replacefield(cell2well2tbl, {'wells', 'wells1'});
    [~, cell2well2tbl] = setupTableMapping(cell2well2tbl, cellwelltbl, {{'cells2', ...
                        'cells'}});
    cell2well2tbl = replacefield(cell2well2tbl, {'wells', 'wells2'});
    map = setupTableMapping(cell2tbl, cell2well2tbl, {'cells1', 'cells2'});
    b = map*b; % b belongs to cell2well2tbl
    
    invb = 1./b;
    
    repfds = {{'cells', 'cells1'}, ...
              {'nodes', 'nodes1'}, ...
              {'wells', 'wells1'}};
    cell_1node_1well_1tbl = replacefield(cellnodewelltbl, repfds);
    repfds = {{'cells', 'cells2'}, ...
              {'nodes', 'nodes2'}, ...
              {'wells', 'wells2'}};
    cell_2node_2well_2tbl = replacefield(cellnodewelltbl, repfds);
    
    [~, prodtbl] = setupTableMapping(cell2well2tbl, cell_1node_1well_1tbl, ...
                                                   {'cells1', 'wells1'});
    [~, prodtbl] = setupTableMapping(prodtbl, cell_2node_2well_2tbl, ...
                                                   {'cells2', 'wells2'});
    map1 = setupTableMapping(cell2well2tbl, prodtbl, {'cells1', 'cells2', ...
                        'wells1', 'wells2'});
    map2 = setupTableMapping(cell_1node_1well_1tbl, prodtbl, {'cells1', 'nodes1', ...
                        'wells1'});    
    map3 = setupTableMapping(cell_2node_2well_2tbl, prodtbl, {'cells2', 'nodes2', ...
                        'wells2'});
    d = (map1*invb).*(map2*a).*(map3*a);
    
    prod_node2tbl = projTable(prodtbl, {'nodes1', 'nodes2'});
    map = setupTableMapping(prodtbl, prod_node2tbl, {'nodes1', 'nodes2'});
    d = map*d; % d belongs to prod_node2tbl

    well2tbl = projTable(cell2well2tbl, {'wells1', 'wells2'});
    map = setupTableMapping(cell2well2tbl, well2tbl, {'wells1', 'wells2'});
    e = map*b;
    

    tbl = prod_node2tbl; % alias
    Ann = sparse(tbl.nodes1, tbl.nodes2, d, nn, nn);
    tbl = nodewelltbl; % alias    
    Anw = sparse(tbl.nodes, tbl.wells, c, nn, nw);
    Awn = Anw';
    tbl = well2tbl; %alias
    Aww = sparse(tbl.wells1, tbl.wells2, e, nw, nw);
    
    Awq = speye(nw);
    Anq = sparse(nn, nw);
    Aext = [[A + Ann, Anw,  Anq];
            [    Awn, Aww, -Awq];
           ];
    
    % Setup control
    ctypetbl.ctypes = (1 : 2)';
    ctypetbl.num    = 2;
    [~, allwellctypetbl] = setupTableMapping(ctypetbl, welltbl, []);
    allwellctypetbl = sortTable(allwellctypetbl, {'ctypes', 'wells'});
    allwellctypetbl = addLocInd(allwellctypetbl, 'wcind');
    
    wellctypetbl.wells = (1 : nw)';
    wellctypetbl.ctypes = welltypes;
    wellctypetbl.num = nw;
    [~, wellctypetbl] = setupTableMapping(allwellctypetbl, wellctypetbl, ...
                                                        {'wells', 'ctypes'});

    tbl = wellctypetbl; %alias
    Acw = sparse(tbl.wells, tbl.wcind, 1, nw, allwellctypetbl.num);

    Aext = [Aext;
           [zeros(nw, nn), Acw]];

    rhs = [zeros(nn + nw, 1); ...
           wellvals];
    
    x = Aext\rhs;
    
    pn = x(1 : nn);
    
end


function c = getWellControlType(W)
    switch W.type
      case 'rate'
        c = 2;
      case 'bhp'
        c = 1;
      otherwise
        error('well type not recognized.');
    end
end
