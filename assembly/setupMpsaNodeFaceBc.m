function [D, bcvals] = setupMpsaNodeFaceBc(bc, G, nnodesperface, tbls)
% The structure bc gives conditions on the nodeface displacement

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


    assert(isfield(bc, 'bcnodefacetbl'), ['this function is meant to set ' ...
                        'boundary conditions at the nodeface'])
    
    nodefacetbl    = tbls.nodefacetbl;
    nodefacecoltbl = tbls.nodefacecoltbl;
    coltbl         = tbls.coltbl;
    facetbl        = tbls.facetbl;
    
    bcnodefacetbl = bc.bcnodefacetbl;
    bcnodefacetbl = bcnodefacetbl.addLocInd('bcinds');
    bcnodefacecoltbl = crossIndexArray(bcnodefacetbl, coltbl, {}, 'optpureproduct', true);
    
    linform = bc.linform;
    linform = reshape(linform', [], 1);
    % linform belongs to bcnodefacecoltbl    

    % We compute the (pseudo) area of the nodeface
    faces = facetbl.get('faces');
    nfareas = 1./nnodesperface.*(G.faces.areas(faces));
    
    % We weight the linear form with the area
    prod = TensorProd();
    prod.tbl1 = facetbl;
    prod.tbl2 = bcnodefacecoltbl;
    prod.tbl3 = bcnodefacecoltbl;
    prod.mergefds = {'faces'};
    prod = prod.setup();
    
    linform = prod.eval(nfareas, linform);
    
    prod = TensorProd();
    prod.tbl1 = facetbl;
    prod.tbl2 = bcnodefacetbl;
    prod.tbl3 = bcnodefacetbl;
    prod.mergefds = {'faces'};
    prod = prod.setup();
    
    bcvals = bc.linformvals;
    bcvals = prod.eval(nfareas, bcvals);
    
    prod = TensorProd();
    prod.tbl1 = bcnodefacecoltbl;
    prod.tbl2 = bcnodefacecoltbl;
    prod.tbl3 = bcnodefacetbl;
    prod.mergefds = {'nodes', 'faces', 'bcinds'};
    prod.reducefds = {'coldim'};
    prod = prod.setup();
     
    D_T = SparseTensor('matlabsparse', true);
    D_T = D_T.setFromTensorProd(linform, prod);

    map = TensorMap();
    map.fromTbl = nodefacecoltbl;
    map.toTbl = bcnodefacecoltbl;
    map.mergefds = {'nodes', 'faces', 'coldim'};
    map = map.setup();
    
    M_T = SparseTensor('matlabsparse', true);
    M_T = M_T.setFromTensorMap(map);
    
    D_T = D_T*M_T;
    
    D = D_T.getMatrix();
    D = D';
    
end


