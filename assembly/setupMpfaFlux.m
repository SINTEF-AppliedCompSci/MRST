function mpfaKgrad = setupMpfaFlux(G, assembly, tbls)
% we set up the mappings
%
%   F1 : nodefacetbl -> intfacetbl 
%   F2 : celltbl -> intfacetbl
%
%  such that the flux u in intfacetbl (interior faces) is given by
%
%  u  = [F1  F2] * [ pnf (pressure at nodefacetbl);
%                    pc  (pressure at celltbl)];
%
%  Then, we proceed with the reduction to remove dependency in pnf
%  (pressure at nodefacetbl)

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

    
    matrices = assembly.matrices;
    nKg = assembly.nKg;
    
    invA11 = matrices.invA11;
    A12    = matrices.A12;
    
    celltbl      = tbls.celltbl;
    cellfacetbl      = tbls.cellfacetbl;
    facetbl          = tbls.facetbl;
    cellnodeface2tbl = tbls.cellnodeface2tbl;
    nodefacetbl      = tbls.nodefacetbl;

    intfaces = find(all(G.faces.neighbors, 2));
    intfacetbl.faces = intfaces;
    intfacetbl = IndexArray(intfacetbl);
    
    map = TensorMap();
    map.fromTbl  = cellfacetbl;
    map.toTbl    = facetbl;
    map.mergefds = {'faces'};
    map = map.setup();
    
    ncellperface = map.eval(ones(cellfacetbl.num, 1));
    
    cno = cellfacetbl.get('cells');
    fno = cellfacetbl.get('faces');
    sgn = 2*(cno == G.faces.neighbors(fno, 1)) - 1;

    prod = TensorProd();
    prod.tbl1 = facetbl;
    prod.tbl2 = cellfacetbl;
    prod.tbl3 = cellfacetbl;
    prod.mergefds = {'faces'};
    prod = prod.setup();
    
    wsgn = prod.eval(1./ncellperface, sgn);
    
    prod = TensorProd();
    prod.tbl1 = cellfacetbl;
    prod.tbl2 = cellnodeface2tbl;
    prod.tbl3 = cellnodeface2tbl;
    prod.replacefds1 = {{'faces', 'faces1'}};
    prod.mergefds = {'cells', 'faces1'};
    prod = prod.setup();
    
    wnKg = prod.eval(wsgn, nKg);

    gen = CrossIndexArrayGenerator();
    gen.tbl1 = cellnodeface2tbl;
    gen.tbl2 = intfacetbl;
    gen.replacefds2 = {{'faces', 'faces1'}};
    gen.mergefds = {'faces1'};
    
    cellnodeintface2tbl = gen.eval();
 
    map = TensorMap();
    map.fromTbl = cellnodeface2tbl;
    map.toTbl = cellnodeintface2tbl;
    map.mergefds = {'cells', 'nodes', 'faces1', 'faces2'};
    map = map.setup();
    
    wnKg = map.eval(wnKg);
    
    prod = TensorProd();
    prod.tbl1 = cellnodeintface2tbl;
    prod.tbl2 = nodefacetbl;
    prod.tbl3 = intfacetbl;
    prod.replacefds1 = {{'faces1', 'faces'}};
    prod.replacefds2 = {{'faces', 'faces2'}};
    prod.reducefds = {'nodes', 'faces2'};
    prod = prod.setup();
    
    F1_T = SparseTensor();
    F1_T = F1_T.setFromTensorProd(wnKg, prod);
    
    F1 = F1_T.getMatrix();

    prod = TensorProd();
    prod.tbl1 = cellnodeintface2tbl;
    prod.tbl2 = celltbl;
    prod.tbl3 = intfacetbl;
    prod.replacefds1 = {{'faces1', 'faces'}};
    prod.reducefds = {'cells'};
    prod = prod.setup();
    
    F2_T = SparseTensor();
    % note the minus sign
    F2_T = F2_T.setFromTensorProd(-wnKg, prod);
    
    F2 = F2_T.getMatrix();
    
    % We set up flux operator
    %
    %  F : celltbl -> intfacetbl
    %
    % We assume Neumann boundary condition for flow so that we have
    %
    % [A11, A12] * [ pnf (pressure at nodefacetbl);
    %                pc  (pressure at celltbl)     ]   =   0;  
    %
    
    mpfaKgrad = F2 - F1*invA11*A12;
end
