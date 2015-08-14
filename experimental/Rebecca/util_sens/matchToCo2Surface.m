function val = matchToCo2Surface(sG,surface,G,fluid)
    h = G.cells.H.*sG/(1-fluid.res_water)-surface.h;
    val=sum(h.^2.*G.cells.volumes);
end