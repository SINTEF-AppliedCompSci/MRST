function bG = getbG_BO(model, p, rv, isLiquid)
    if model.vapoil
        bG = model.fluid.bG(p, rv, isLiquid);
    else
        bG = model.fluid.bG(p);
    end
end