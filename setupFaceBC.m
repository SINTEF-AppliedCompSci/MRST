function nodefacebc = setupFaceBC(bc, G, tbls)
    
    nodefacetbl    = tbls.nodefacetbl;
    nodefacecoltbl = tbls.nodefacecoltbl;
    coltbl         = tbls.coltbl;

    % note that bcfacetbl, as constructed below, has repeated indices (same face can
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


