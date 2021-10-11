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
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
    tbls = vagstruct.tbls;
    
    celltbl = tbls.celltbl;
    nodetbl = tbls.nodetbl;    
    cellnodetbl = tbls.cellnodetbl;
    cellnode2tbl = tbls.cellnode2tbl;
    
    nc = celltbl.num;
    nn = nodetbl.num;
    
    %% Assembly of A11 (nodes - nodes)
    node2tbl = projIndexArray(cellnode2tbl, {'nodes1', 'nodes2'});
    
    map = TensorMap();
    map.fromTbl = cellnode2tbl;
    map.toTbl = node2tbl;
    map.mergefds = {'nodes1', 'nodes2'};
    map = map.setup();
    
    A11 = map.eval(Atrans);
    
    % setup matrix from A11
    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = node2tbl;
    map.replaceFromTblfds = {{'nodes', 'nodes1'}};
    map.mergefds = {'nodes1'};
    ind1 = getDispatchInd(map);
    
    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = node2tbl;
    map.replaceFromTblfds = {{'nodes', 'nodes2'}};
    map.mergefds = {'nodes2'};
    ind2 = getDispatchInd(map);

    A11m = sparse(ind1, ind2, A11, nn, nn);
    

    %% Assembly of A12 (nodes - cells)
    
    map = TensorMap();
    map.fromTbl = cellnode2tbl;
    map.toTbl = cellnodetbl;
    map.replaceFromTblfds = {{'nodes1', 'nodes'}};
    map.mergefds = {'cells', 'nodes'};
    map = map.setup();
    
    A12 = -map.eval(Atrans);
    
    % setup matrix from A12
    map = TensorMap();
    map.fromTbl = nodetbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'nodes'};
    ind1 = getDispatchInd(map);
    
    map = TensorMap();
    map.fromTbl = celltbl;
    map.toTbl = cellnodetbl;
    map.mergefds = {'cells'};
    ind2 = getDispatchInd(map);

    A12m = sparse(ind1, ind2, A12, nn, nc);
    

    %% Assembly of A22

    map = TensorMap();
    map.fromTbl = cellnode2tbl;
    map.toTbl = celltbl;
    map.mergefds = {'cells'};
    map = map.setup();
    
    A22 = map.eval(Atrans);
    
    invA22 = 1./A22;
    ind = (1 : nc)';
    invA22m = sparse(ind, ind, invA22, nc, nc);

    % We have (written in matrix form)
    % [[A11 , A12]   *    [[pn]     = [[f]
    %  [A12', A22]]        [pc]]       [g]]
    % 
    % A22 is diagonal and can be inverted directly. We get
    % pc = invA22*g - invA22*A12'*pn
    % and
    % A * pn = rhs
    % for
    % A = A11 - A12*invA22*A12'
    % and
    % rhs = f - A12*invA22*g
    %
    
    % system matrix:
    A = A11m - A12m*invA22m*A12m';

    rhsfun = @(f, g) assembleRHS(f, g, A12m, invA22m);
    computeCellPressure = @(pn, g) assembleCellPressure(pn, g, A12m, invA22m);
   
    matrices = struct('A11', A11, ...
                      'A12', A12, ...
                      'A22', A22);

    tbls.node2tbl = node2tbl;
    
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
