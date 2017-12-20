function frac = getNonUnitMassFraction(model, molfraction, molfraction_normalized)
% Internal utility. Intentionally undocumented.
    model = model.EOSModel;
    ncomp = numel(molfraction);
    mass = cell(1, ncomp);
    total = 0;
    mw = model.fluid.molarMass;
    if iscell(molfraction)
        for i = 1:numel(molfraction)
            mass{i} = molfraction{i}.*mw(i);
            total = total + molfraction_normalized{i}.*mw(i);
        end
        frac = cellfun(@(x) x./total, mass, 'UniformOutput', false);
    else
        mass = molfraction.*mw;
        total = sum(molfraction_normalized.*mw, 2);
        frac = mass./total;
    end
    
end