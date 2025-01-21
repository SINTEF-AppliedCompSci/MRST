function [nodefacebc, tbls, mappings] = setupFaceBC2(bc, G, tbls, mappings, varargin)

    opt = struct('useVirtual', false);
    opt = merge_options(opt, varargin{:});
    
    useVirtual = opt.useVirtual;
    
    nodefacetbl    = tbls.nodefacetbl;
    nodefacevectbl = tbls.nodefacevectbl;
    vectbl         = tbls.vectbl;

    % Note that bcfacetbl, as constructed below, has repeated indices (same face can
    % have several linear forms which are imposed on it). We add a local index
    % denoted bcinds, which makes all those (now multiple-)index unique.
    bcfacetbl.faces = bc.extfaces;
    bcfacetbl = IndexArray(bcfacetbl);
    
    bcfacetbl = bcfacetbl.addLocInd('bcinds');
    
    [bcnodefacetbl, indstruct] = crossIndexArray(bcfacetbl, nodefacetbl, {'faces'});
    bcface_from_bcnodeface   = indstruct{1}.inds;
    nodeface_from_bcnodeface = indstruct{2}.inds;
    
    bcfacevectbl     = crossIndexArray(bcfacetbl    , vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);
    bcnodefacevectbl = crossIndexArray(bcnodefacetbl, vectbl, {}, 'optpureproduct', true, 'virtual', useVirtual);
    
    bcvals = bc.linformvals;
    % bcvals belongs to bcfacetbl;
    map = TensorMap();
    map.fromTbl  = bcfacetbl;
    map.toTbl    = bcnodefacetbl;
    map.mergefds = {'faces', 'bcinds'};

    if useVirtual

        map.pivottbl = bcnodefacetbl;
        map.dispind1 = bcface_from_bcnodeface;
        map.dispind2 = (1 : bcnodefacetbl.num)';
        map.issetup = true;
        
    else
    
        map = map.setup();
    end
    
    bcvals = map.eval(bcvals);
    % bcvals belongs to bcnodefacetbl;
    
    
    linform = bc.linform;
    linform = reshape(linform', [], 1);
    % linform belongs to bcfacevectbl    
    
    map = TensorMap();
    map.fromTbl  = bcfacevectbl;
    map.toTbl    = bcnodefacevectbl;
    map.mergefds = {'faces', 'bcinds', 'vec'};

    if useVirtual

        map.pivottbl = bcnodefacevectbl;

        N = bcnodefacetbl.num;
        [vec, i] = ind2sub([vectbl.num], (1 : bcnodefacevectbl.num)');
        map.dispind1 = sub2ind([vectbl.num,  N], vec, bcface_from_bcnodeface(i));
        map.dispind2 = (1 : bcnodefacevectbl.num)';
        map.issetup = true;
        
    else
    
        map = map.setup();
    end
    
    linform = map.eval(linform);
    % linform now belongs to bcnodefacevectbl
    
    bcnodefacetbl = replacefield(bcnodefacetbl, {{'bcinds', ''}});
    nodefacebc.bcnodefacetbl = bcnodefacetbl;
    nodefacebc.linform       = linform;
    nodefacebc.linformvals   = bcvals;

    tbls.bcnodefacetbl = bcnodefacetbl;
    
    mappings.bcface_from_bcnodeface   = bcface_from_bcnodeface;
    mappings.nodeface_from_bcnodeface = nodeface_from_bcnodeface;
    
end

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

