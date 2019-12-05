function [A, operators] = setupSystem(vagstruct, G)
%
%
% SYNOPSIS:
%   function [A, operators] = setupSystem(vagstruct, G)
%
% DESCRIPTION:
%
% PARAMETERS:
%   vagstruct - as returned by computeVagTrans
%   G         - Grid
%
% RETURNS:
%   A         - System matrix with nodal degree of freedom
%   operators - structure containing
%
%               operators.rhsfun   - function to compute the right-hand side:
%                                    The function takes two vector flux values, corresponding to
%                                    flux at nodes and cells, see function `assembleRHS` at the end of this file 
%
%               operators.matrices - intermediate matrices that are assembled
%
%               operators.tbls     - tables describing various connectivities
%
%               operators.computeCellPressure - function that computes the cell pressures given the nodal pressures.
%
% EXAMPLE:
%   `linearpressuretest`
%
% SEE ALSO:
%   `computeVagTrans`, `incompVAG`

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
    
    Atrans = vagstruct.A;
    cellnode2tbl = vagstruct.cellnode2tbl;
    
    nc = G.cells.num;
    clear celltbl;
    celltbl.cells = (1 : nc)';
    celltbl.num   = nc;
    celltbl       = addLocInd(celltbl, 'cind');

    nn = G.nodes.num;
    clear nodetbl;
    nodetbl.nodes = (1 : nn)';
    nodetbl.num   = nn;
    nodetbl       = addLocInd(nodetbl, 'nind');

    %% Assembly of A11
    node2tbl = projTable(cellnode2tbl, {'nodes1', 'nodes2'});
    [~, node2tbl] = setupTableMapping(node2tbl, nodetbl, {{'nodes1', 'nodes'}});
    node2tbl = replacefield(node2tbl, {'nind', 'nind1'});
    [~, node2tbl] = setupTableMapping(node2tbl, nodetbl, {{'nodes2', 'nodes'}});
    node2tbl = replacefield(node2tbl, {'nind', 'nind2'});

    map = setupTableMapping(cellnode2tbl, node2tbl, {'nodes1', 'nodes2'});
    A11 = map*Atrans;
    tbl = node2tbl;
    A11m = sparse(tbl.nind1, tbl.nind2, A11, nodetbl.num, nodetbl.num);

    %% Assembly of A12
    cellnodetbl = projTable(cellnode2tbl, {'cells', 'nodes1'});
    [~, cellnodetbl] = setupTableMapping(cellnodetbl, nodetbl, {{'nodes1', 'nodes'}});
    cellnodetbl = replacefield(cellnodetbl, {'nind', 'nind1'});

    map = setupTableMapping(cellnode2tbl, cellnodetbl, {'cells', 'nodes1'});
    A12 = -map*Atrans;
    tbl = cellnodetbl;
    A12m = sparse(tbl.nind1, tbl.cells, A12, nodetbl.num, celltbl.num);

    %% Assembly of A12

    cell2tbl = duplicatefield(celltbl, {'cells', {'cells1', 'cells2'}});
    map = setupTableMapping(cellnode2tbl, cell2tbl, {{'cells', 'cells1'}});
    A22 = map*Atrans;
    invA22 = 1./A22;
    tbl = cell2tbl;
    invA22m = sparse(tbl.cells1, tbl.cells2, invA22, tbl.num, tbl.num);

    % We have
    %
    % [[A11 , A12]  * [[pn]  = [[f]
    %  [A12', A22]]    [pc]]    [g]] 
    % 
    % A22 is diagonal and can be inverted directly. We get 
    %
    % pc = invA22*g - invA22*A12'*pn
    %
    % and
    %
    % A * pn = rhs
    %
    % for 
    %
    % A = A11 - A12*invA22*A12'
    %
    % and 
    %
    % rhs = f - A12*invA22*g
    %
    
    % system matrix:
    A = A11m - A12m*invA22m*A12m';

    rhsfun = @(f, g) assembleRHS(f, g, A12m, invA22m);
    computeCellPressure = @(pn, g) assembleCellPressure(pn, g, A12m, invA22m);
   
    matrices = struct('A11', A11, ...
                      'A12', A12, ...
                      'A22', A22);

    % cleaning the tables before exporting
    cellnodetbl = rmfield(cellnodetbl, 'nind1');
    cellnodetbl = replacefield(cellnodetbl, {'nodes1', 'nodes'});
    node2tbl = rmfield(node2tbl, 'nind1');
    node2tbl = rmfield(node2tbl, 'nind2');
    
    tbls = struct('node2tbl'   , node2tbl   , ...
                  'cellnodetbl', cellnodetbl, ...
                  'cell2tbl'   , cell2tbl);
    
    operators.rhsfun   = rhsfun;
    operators.matrices = matrices;
    operators.tbls     = tbls;
    operators.computeCellPressure = computeCellPressure;
   
end


function rhs = assembleRHS(f, g, A12, invA22)
% f :  flux source term at nodes
% g :  flux source term at cells
% A12, invA2 : Assembly matrices as computed above
        rhs = f - A12*invA22*g;
end

function pc = assembleCellPressure(pn, g, A12, invA22)
% pn : nodal pressure
% g :  flux source term at cells
% A12, invA2 : Assembly matrices as computed above
    pc = invA22*g - invA22*A12'*pn;
end
