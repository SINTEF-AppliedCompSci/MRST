function c = getEOSComponent(model, p, T, name, ci)
    names_hc = model.EOSModel.fluid.names;
    hcpos = strcmp(names_hc, name);
    n_hc = numel(names_hc);
    z = zeros(1, n_hc);
    z(hcpos) = 1;
    [L, ~, ~, ~, ~, rhoL, rhoV] = standaloneFlash(p, T, z, model.EOSModel);
    Lm = L.*rhoL./(rhoL.*L + rhoV.*(1-L));
    Vm = 1 - Lm;
    if model.water
        frac = [0, Lm, Vm];
        rho = [model.fluid.rhoWS, rhoL, rhoV];
    else
        frac = [Lm, Vm];
        rho = [rhoL, rhoV];
    end
    c = EquationOfStateComponent(name, p, T, ci, frac, rho, model.EOSModel.fluid.molarMass(hcpos));
end