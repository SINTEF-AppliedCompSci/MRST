function frac = getNonUnitMoleFraction(model, massfraction, massfraction_normalized)
% Internal utility. Intentionally undocumented.
    model = model.EOSModel;
    ncomp = numel(massfraction);
    mass = cell(1, ncomp);
    total = 0;
    mw = model.fluid.molarMass;
    if iscell(massfraction)
        for i = 1:numel(massfraction)
            mass{i} = massfraction{i}./mw(i);
            total = total + massfraction_normalized{i}./mw(i);
        end
        frac = cellfun(@(x) x./total, mass, 'UniformOutput', false);
    else
        mass = massfraction./mw;
        total = sum(massfraction_normalized./mw, 2);
        frac = mass./total;
    end
    
end