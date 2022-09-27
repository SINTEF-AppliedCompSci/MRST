function assembly = assembleMPSA(G, prop, loadstruct, eta, tbls, mappings, varargin)
% Assembly of MPSA-weak
%
% Reference paper:
% Finite volume methods for elasticity with weak symmetry
% Keilegavlen, Eirik and Nordbotten, Jan Martin
% International Journal for Numerical Methods in Engineering
% 2017

    % the solution is given by the system
    %
    % A = [[A11, A12, -D];
    %      [A21, A22,  0];
    %      [D' , 0  ,  0]];
    %
    % u = [u (displacement at nodefacecoltbl);
    %      u (displacement at cellcoltbl);
    %      lagmult (forces in the linear directions at the boundary)];
    %
    % f = [extforce (force at nodefacecoltbl);
    %      force    (volumetric force at cellcoltbl);
    %      bcvals   (for the linear form at the boundary)];
    %
    % A*u = f
    %
    % Note: extforce is sparse and should only give contribution at facets
    % that are at the boundary
    %
    % By construction of the method, the matrix A11 is block-diagonal. Hence,
    % we invert it directly and reduce to a cell-centered scheme.

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

    
    opt = struct('bcetazero'       , true , ...
                 'assemblyMatrices', false, ...
                 'addAdOperators'  , false, ...
                 'extraoutput'     , false);
    opt = merge_options(opt, varargin{:});
    
    % Recover IndexArrays
    coltbl                = tbls.coltbl;
    celltbl               = tbls.celltbl;
    facetbl               = tbls.facetbl;
    nodetbl               = tbls.nodetbl;
    cellnodetbl           = tbls.cellnodetbl;
    nodefacetbl           = tbls.nodefacetbl;
    cellcoltbl            = tbls.cellcoltbl;
    nodecoltbl            = tbls.nodecoltbl;
    nodefacecoltbl        = tbls.nodefacecoltbl;
    cellnodefacetbl       = tbls.cellnodefacetbl;
    cellnodecoltbl        = tbls.cellnodecoltbl;
    cellnodecolrowtbl     = tbls.cellnodecolrowtbl;
    cellnodefacecoltbl    = tbls.cellnodefacecoltbl;
    cellnodefacecolrowtbl = tbls.cellnodefacecolrowtbl;
    colrowtbl             = tbls.colrowtbl;
    nodecolrowtbl         = tbls.nodecolrowtbl;
    col2row2tbl           = tbls.col2row2tbl;
    cellcol2row2tbl       = tbls.cellcol2row2tbl;
    cellnodecol2row2tbl   = tbls.cellnodecol2row2tbl;
    
    cell_from_cellnode         = mappings.cell_from_cellnode;
    node_from_cellnode         = mappings.node_from_cellnode;
    cellnode_from_cellnodeface = mappings.cellnode_from_cellnodeface;
    nodeface_from_cellnodeface = mappings.nodeface_from_cellnodeface;
    
    % Some shortcuts
    c_num     = celltbl.num;
    n_num     = nodetbl.num;
    cnf_num   = cellnodefacetbl.num;
    cnfc_num  = cellnodefacecoltbl.num;
    cn_num    = cellnodetbl.num;
    cncr_num  = cellnodecolrowtbl.num;
    nf_num    = nodefacetbl.num;
    nfc_num   = nodefacecoltbl.num;
    cnfcr_num = cellnodefacecolrowtbl.num;
    d_num     = coltbl.num;
    
    dim = coltbl.num;

    C = setupStiffnessTensor(prop, tbls);
    
    map = TensorMap();
    map.fromTbl = nodefacetbl;
    map.toTbl = facetbl;
    map.mergefds = {'faces'};
    map = map.setup();
    
    nnodesperface = map.eval(ones(nodefacetbl.num, 1));
    
    bc = loadstruct.bc;
    opts = struct('eta', eta, ...
                  'bcetazero', opt.bcetazero);
    [matrices, bcvals, extra] = coreMpsaAssembly(G, C, bc, nnodesperface, tbls, mappings, opts);
    
    extforce = loadstruct.extforce;
    force = loadstruct.force;

    fullrhs{1} = extforce;
    fullrhs{2} = force;
    fullrhs{3} = bcvals;
    
    matrices.fullrhs = fullrhs;
    
    % We reduced the system (shur complement) using invA11
    % We obtain system of the form
    %
    % B*u = rhs
    %
    % where
    %
    % B = [[B11, B12];
    %      [B21, B22]];
    %
    % u = [u (displacement at cellcoltbl);
    %      lagmult];
    %
    % rhs = [-A21*invA11*extforce;  +  [force;
    %        -D'*invA11*extforce  ]     bcvals]
    
    invA11 = matrices.invA11; 
    A12 = matrices.A12;
    A21 = matrices.A21;
    A22 = matrices.A22;
    D = matrices.D;
    
    B11 = A22 - A21*invA11*A12;
    B12 = A21*invA11*D;
    B21 = -D'*invA11*A12;
    B22 = D'*invA11*D;
    
    B = [[B11, B12]; ...
         [B21, B22]];
    
    adrhs{1} = -A21*invA11*extforce + force; 
    adrhs{2} = -D'*invA11*extforce + bcvals;
    
    rhs = vertcat(adrhs{:});

    % Assembly of operator to compute u_{nodefacecoltbl} from solution of the system
    % (which consists of the concatenation of u_{cellcol} and lagmult) and
    % extforce which is a force in nodefacecoltbl
    %
    % We have  u_{nodefacecoltbl} = R1*sol + R2*extforce
    %
    % where R1 = invA11*[-A12, D] and R2 = invA11
    
    R1 = invA11*[-A12, D];
    R2 = invA11;
    g = extra.g;
    
    assembly = struct('B'       , B       , ...
                      'rhs'     , rhs     , ...
                      'g'       , g       , ...
                      'extforce', extforce, ...
                      'R1'      , R1      , ...
                      'R2'      , R2);
    
    if opt.assemblyMatrices
        assembly.matrices = matrices;
    end
    
    if opt.addAdOperators

        adB = cell(2, 2);
        adB{1, 1} = B11;
        adB{2, 1} = B21;
        adB{1, 2} = B12;
        adB{2, 2} = B22;
        
        % Setup divergence operator (dilatation)
        div = matrices.div;
        divop = @(sol) mpsaDivOperator(sol, extforce, R1, R2, div);

        % Set up consistent divergence operator (not efficient implementation)
        alpha = ones(G.cells.num, 1);
        % we need to send number of nodes per cell
        map = TensorMap();
        map.fromTbl = cellnodetbl;
        map.toTbl = celltbl;
        map.mergefds = {'cells'};
        map = map.setup();
        nnodespercell = map.eval(ones(cellnodetbl.num, 1));
        cassembly = assembleCouplingTerms(G, eta, alpha, nnodespercell, tbls, mappings);
        cdiv{1} = cassembly.divconsnf;
        cdiv{2} = cassembly.divconsc;
        cdivop = @(unf, uc) cdivopFunc(unf, uc, cdiv);
        
        % Setup node face displacement operator
        fndisp{1} = -invA11*A12;
        fndisp{2} = invA11*D;
        fndisp{3} = invA11;
        fndispop = @(u, lm) facenodedispopFunc(u, lm, extforce, fndisp);
        
        % Setup stress operator
        aver = cellAverageOperator(tbls, mappings);
        stress{1} = matrices.C1;
        stress{2} = matrices.C2;
        stressop = @(unf, uc) stressopFunc(unf, uc, stress, aver);
        
        adoperators.B     = adB;
        adoperators.rhs   = adrhs;        
        
        adoperators.divop          = divop;
        adoperators.cdivop         = cdivop;
        adoperators.facenodedispop = fndispop;
        adoperators.stressop       = stressop;
        
        assembly.adoperators = adoperators;
        
    end
        
    if opt.extraoutput
        div = matrices.div;
        assembly.divop = @(sol) mpsaDivOperator(sol, extforce, R1, R2, div);
    end
    
end

function fndisp = facenodedispopFunc(u, lm, extforce, fndisp)
    fndisp = fndisp{1}*u + fndisp{2}*lm + fndisp{3}*extforce;
end

function stress = stressopFunc(unf, uc, stress, aver)
    % Stress at each cell-node region (corner)
    stress = stress{1}*unf + stress{2}*uc;
    % Compute cell average
    stress = aver*stress;
end

function cdiv = cdivopFunc(unf, uc, cdiv)
    cdiv = cdiv{1}*unf + cdiv{2}*uc;
end

