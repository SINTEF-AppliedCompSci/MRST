function output = setupMpsaNodeFaceBc(bc, G, nnodesperface, tbls, mappings, varargin)
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

    opt = struct('useVirtual', false);
    opt = merge_options(opt, varargin{:});

    useVirtual = opt.useVirtual;
    
    assert(isfield(bc, 'bcnodefacetbl'), ['this function is meant to set ' ...
                        'boundary conditions at the nodeface'])
    
    nodefacetbl    = tbls.nodefacetbl;
    nodefacevectbl = tbls.nodefacevectbl;
    vectbl         = tbls.vectbl;
    facetbl        = tbls.facetbl;
    
    bcnodefacetbl = bc.bcnodefacetbl;

    map = TensorMap();
    map.fromTbl  = nodefacetbl;
    map.toTbl    = bcnodefacetbl;
    map.mergefds = {'nodes', 'faces'};
    nodeface_from_bcnodeface = map.getDispatchInd();

    bcnodefacetbl = bcnodefacetbl.addLocInd('bcinds');

    bcnodefacevectbl = crossIndexArray(bcnodefacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);
    
    linform = bc.linform;
    linform = reshape(linform', [], 1);
    % linform belongs to bcnodefacevectbl    

    % We compute the (pseudo) area of the nodeface
    faces = facetbl.get('faces');
    nfareas = 1./nnodesperface.*(G.faces.areas(faces));
    
    % We weight the linear form with the area
    prod = TensorProd();
    prod.tbl1 = facetbl;
    prod.tbl2 = bcnodefacevectbl;
    prod.tbl3 = bcnodefacevectbl;
    prod.mergefds = {'faces'};

    if useVirtual
        
        prod.pivottbl = bcnodefacevectbl;
        
        [vec, i] = ind2sub([vectbl.num, bcnodefacetbl.num], (1 : bcnodefacevectbl.num)');

        prod.dispind1 = mappings.face_from_nodeface(nodeface_from_bcnodeface(i));
        prod.dispind2 = (1 : bcnodefacevectbl.num)';
        prod.dispind3 = (1 : bcnodefacevectbl.num)';
        
        prod.issetup = true;

    else
        prod = prod.setup();
    end
    
    linform = prod.eval(nfareas, linform);
    
    prod = TensorProd();
    prod.tbl1 = facetbl;
    prod.tbl2 = bcnodefacetbl;
    prod.tbl3 = bcnodefacetbl;
    prod.mergefds = {'faces'};

    if useVirtual
        
        prod.pivottbl = bcnodefacetbl;

        i = (1 : bcnodefacetbl.num)';
        j = nodeface_from_bcnodeface(i);
        j = mappings.face_from_nodeface(j);
        prod.dispind1 = j;
        prod.dispind2 = i;
        prod.dispind3 = i;
        
        prod.issetup = true;

    else
        prod = prod.setup();
    end
    
    bcvals = bc.linformvals;
    bcvals = prod.eval(nfareas, bcvals);
    
    prod = TensorProd();
    prod.tbl1      = bcnodefacevectbl;
    prod.tbl2      = bcnodefacevectbl;
    prod.tbl3      = bcnodefacetbl;
    prod.mergefds  = {'nodes', 'faces', 'bcinds'};
    prod.reducefds = {'vec'};

    if useVirtual
        
        prod.pivottbl = bcnodefacevectbl;

        [vec, i] = ind2sub([vectbl.num, bcnodefacetbl.num], (1 : bcnodefacevectbl.num)');

        prod.dispind1 = (1 : bcnodefacevectbl.num)';
        prod.dispind2 = (1 : bcnodefacevectbl.num)';
        prod.dispind3 = i;
        
        prod.issetup = true;

    else
        prod = prod.setup();
    end
     
    D = prod.setupMatrix(linform);

    map = TensorMap();
    map.fromTbl  = nodefacevectbl;
    map.toTbl    = bcnodefacevectbl;
    map.mergefds = {'nodes', 'faces', 'vec'};

    if useVirtual
        
        map.pivottbl = bcnodefacevectbl;
        
        [vec, i] = ind2sub([vectbl.num, bcnodefacetbl.num], (1 : bcnodefacevectbl.num)');
        
        map.dispind1 = sub2ind([vectbl.num, nodefacetbl.num], vec, nodeface_from_bcnodeface(i));
        map.dispind2 = (1 : bcnodefacevectbl.num)';
        map.issetup = true;
    else
        map = map.setup();
    end
    
    M = map.getMatrix();
    
    D = D*M;
    
    D = D';

    tbls.bcnodefacetbl    = bcnodefacetbl;
    tbls.bcnodefacevectbl = bcnodefacevectbl;

    mappings.nodeface_from_bcnodeface = nodeface_from_bcnodeface;
    
    output = struct('D'       ,D     , ...
                    'bcvals'  ,bcvals, ...
                    'tbls'    ,tbls  , ...
                    'mappings', mappings);
    
end


