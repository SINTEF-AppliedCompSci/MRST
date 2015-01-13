function bG = getbO_BO(model, p, rs, isSaturated)
    if isprop(model, 'disgas') && model.disgas
        bG = model.fluid.bO(p, rs, isSaturated);
    else
        bG = model.fluid.bO(p);
    end
end