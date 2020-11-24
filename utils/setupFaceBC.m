function nodefacebc = setupFaceBC(bc, G, tbls)
    
    nodefacetbl    = tbls.nodefacetbl;
    nodefacecoltbl = tbls.nodefacecoltbl;
    coltbl         = tbls.coltbl;

    % Note that bcfacetbl, as constructed below, has repeated indices (same face can
    % have several linear forms which are imposed on it). We add a local index
    % denoted bcinds, which makes all those (now multiple-)index unique.
    bcfacetbl.faces = bc.extfaces;
    bcfacetbl = IndexArray(bcfacetbl);
    
    bcfacetbl = bcfacetbl.addLocInd('bcinds');
    
    bcfacecoltbl = crossIndexArray(bcfacetbl, coltbl, {}, 'optpureproduct', true);
    
    bcnodefacetbl = crossIndexArray(bcfacetbl, nodefacetbl, {'faces'});
    bcnodefacecoltbl = crossIndexArray(bcnodefacetbl, coltbl, {}, 'optpureproduct', ...
                                  true);
    
    bcvals = bc.linformvals;
    % bcvals belongs to bcfacetbl;
    map = TensorMap();
    map.fromTbl = bcfacetbl;
    map.toTbl = bcnodefacetbl;
    map.mergefds = {'faces', 'bcinds'};
    map = map.setup();
    
    bcvals = map.eval(bcvals);
    % bcvals belongs to bcnodefacetbl;
    
    
    linform = bc.linform;
    linform = reshape(linform', [], 1);
    % linform belongs to bcfacecoltbl    
    
    map = TensorMap();
    map.fromTbl = bcfacecoltbl;
    map.toTbl = bcnodefacecoltbl;
    map.mergefds = {'faces', 'bcinds', 'coldim'};
    map = map.setup();
    
    linform = map.eval(linform);
    % linform now belongs to bcnodefacecoltbl
    
    bcnodefacetbl = replacefield(bcnodefacetbl, {{'bcinds', ''}});
    nodefacebc.bcnodefacetbl = bcnodefacetbl;
    nodefacebc.linform       = linform;
    nodefacebc.linformvals   = bcvals;

    
end

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
