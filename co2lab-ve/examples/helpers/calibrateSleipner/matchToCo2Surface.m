function val = matchToCo2Surface(sG,surface,G,fluid)
% "surface.h" is permitted to contain nan's (i.e., missing data), in which
% case the following calculation omits nans in order to quantify the match
% between simulated heights and (present) data heights

    h = G.cells.H(~isnan(surface.h)) .* sG(~isnan(surface.h)) / (1-fluid.res_water) ...
        - surface.h(~isnan(surface.h));
    val = sum( h.^2 .* G.cells.volumes(~isnan(surface.h)) );
    
    %h = G.cells.H.*sG/(1-fluid.res_water)-surface.h;
    %val=sum(h.^2.*G.cells.volumes, 'omitnan'); %@@ doesn't work for ADIs
end