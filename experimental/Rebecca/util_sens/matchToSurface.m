function val = matchToCo2Surface(state,surface,G,fluid)
    h = state.sG/(1-fluid.res_water)-surface.h;
    val=sum(h.^2.*Gt.cells.volumes);
end