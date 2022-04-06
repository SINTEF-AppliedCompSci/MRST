function RUNSPEC = generateRUNSPEC(model)
    RUNSPEC = struct();
    if model.oil
        RUNSPEC.OIL = true;
    end
    if model.gas
        RUNSPEC.GAS = true;
    end
    if model.water
        RUNSPEC.WATER = true;
    end
    if isprop(model, 'disgas') && model.disgas
        RUNSPEC.DISGAS = true;
    end
    if isprop(model, 'vapoil') && model.vapoil
        RUNSPEC.VAPOIL = true;
    end
end
