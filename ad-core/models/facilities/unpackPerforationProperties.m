function [p, mob, rho, dissolved, comp, wellvars] = unpackPerforationProperties(packed)
    p           = packed.pressure;
    mob         = packed.mob;
    rho         = packed.rho;
    dissolved   = packed.dissolved;
    comp        = packed.components;
    wellvars    = packed.extravars;
end