function assembly = assembleCouplingTerms(G, eta, alpha, nnodespercell, tbls, mappings)
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


    cellnodefacevectbl = tbls.cellnodefacevectbl;
    nodefacevectbl     = tbls.nodefacevectbl;
    cellnodefacetbl    = tbls.cellnodefacetbl;
    cellnodetbl        = tbls.cellnodetbl;
    cellvectbl         = tbls.cellvectbl;
    celltbl            = tbls.celltbl;
    
    % We fetch the vector g, which belongs to cellnodefacevectbl and is used to
    % construct the consistent divergence operator.
    g = computeConsistentGradient2(G, eta, tbls, mappings);

    % We fetch the vector facetNormals, which belongs to cellnodefacevectbl and is
    % used to construct the finite volume divergence operator.
    normals = computeFacetNormals(G, cellnodefacetbl);

    % Multiply with Biot's coefficient alpha
    prod = TensorProd();
    prod.tbl1 = celltbl;
    prod.tbl2 = cellnodefacevectbl;
    prod.tbl3 = cellnodefacevectbl;
    prod.mergefds = {'cells'};
    prod = prod.setup();
    
    normals = prod.eval(alpha, normals);
    
    % We setup the finite volume divergence operator
    % divfv : nodefacevectbl -> celltbl
    
    prod = TensorProd();
    prod.tbl1 = cellnodefacevectbl;
    prod.tbl2 = nodefacevectbl;
    prod.tbl3 = celltbl;
    prod.reducefds = {'nodes', 'faces', 'vec'};
    prod = prod.setup();
    
    divfv_T = SparseTensor();
    divfv_T = divfv_T.setFromTensorProd(normals, prod);
    divfv = divfv_T.getMatrix();
    
    % We setup the consistent divergence operator
    % It consists of two parts,
    % divconsnf : nodefacevectbl -> celltbl
    % divconsc : celltbl -> celltbl        
    
    % We agregate the contribution at each cell corner.
    % We use equal weights mcoef = (1/(number of nodes per cell)*(volume of the cell)).
    
    mcoef = 1./nnodespercell;    
    
    cno = celltbl.get('cells');
    vols = G.cells.volumes(cno);
    
    mcoef = vols/nnodespercell;
    
    prod = TensorProd();
    prod.tbl1 = celltbl;
    prod.tbl2 = cellnodefacevectbl;
    prod.tbl3 = cellnodefacevectbl;
    prod.mergefds = {'cells'};
    prod = prod.setup();
    
    mg = prod.eval(mcoef, g);
    
    % Multiply with Biot's coefficient alpha
    prod = TensorProd();
    prod.tbl1 = celltbl;
    prod.tbl2 = cellnodefacevectbl;
    prod.tbl3 = cellnodefacevectbl;
    prod.mergefds = {'cells'};
    prod = prod.setup();
    
    mg = prod.eval(alpha, mg);
    
    % We assemble divconsnf : nodefacevectbl -> celltbl
    prod = TensorProd();
    prod.tbl1 = cellnodefacevectbl; 
    prod.tbl2 = nodefacevectbl;
    prod.tbl3 = celltbl;
    prod.reducefds = {'nodes', 'faces', 'vec'};
    prod = prod.setup();
    
    divconsnf_T = SparseTensor();
    divconsnf_T = divconsnf_T.setFromTensorProd(mg, prod);
    divconsnf = divconsnf_T.getMatrix();
    
    % We assemble divconsc : cellvectbl -> celltbl    
    prod = TensorProd();
    prod.tbl1 = cellnodefacevectbl; 
    prod.tbl2 = cellvectbl;
    prod.tbl3 = celltbl;
    prod.mergefds = {'cells'};
    prod.reducefds = {'vec'};
    prod = prod.setup();
    
    divconsc_T = SparseTensor();
    % Beware of minus sign below
    divconsc_T = divconsc_T.setFromTensorProd(-mg, prod);
    divconsc = divconsc_T.getMatrix();
    
    assembly = struct('divfv'    , divfv    , ...
                      'divconsnf', divconsnf, ...
                      'divconsc' , divconsc);
    
end

