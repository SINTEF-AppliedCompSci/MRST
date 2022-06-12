function RUNSPEC = generateRUNSPEC(model, RUNSPEC)
    % Generate RUNSPEC from model. Very limited internal function.
    if nargin == 1
        RUNSPEC = struct();
    end
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
