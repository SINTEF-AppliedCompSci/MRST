function model = convertBlackOilModelToCompositionalModel(bomodel, varargin)
% Convert blackoil model to compositional-like model
    props = convertRSRVtoKValue(bomodel, varargin{:});

    K_o_interp = getMultiDimInterpolator({props.k_values.z_o(end:-1:1, 1), props.k_values.pressure(1, :)'}, props.k_values.K_o(end:-1:1, :));
    K_g_interp = getMultiDimInterpolator({props.k_values.z_g(:, 1), props.k_values.pressure(1, :)'}, props.k_values.K_g);

    K_o = @(p, T, z) K_o_interp(z{1}, p);
    K_g = @(p, T, z) K_g_interp(z{2}, p);

    isSat_interp = getMultiDimInterpolator({props.k_values.z_g(:, 1), props.k_values.pressure(1, :)'}, (props.k_values.freeGas));
    isStableFn = @(p, T, z) isSat_interp(z{2}, p) == 0;
    % isSatFn = @(p, T, z) ~(z{1} > 1e-8 & z{2} > 1e-8);




    muO_interp = getMultiDimInterpolator({props.properties.x_o(:, 1), props.properties.pressure(1, :)'}, 1./props.properties.muO);
    rhoO_interp = getMultiDimInterpolator({props.properties.x_o(:, 1), props.properties.pressure(1, :)'}, props.properties.rhoO);

    muG_interp = getMultiDimInterpolator({props.properties.y_g(:, 1), props.properties.pressure(1, :)'}, 1./props.properties.muG);
    rhoG_interp = getMultiDimInterpolator({props.properties.y_g(:, 1), props.properties.pressure(1, :)'}, props.properties.rhoG);

    muO = @(p, T, x) 1./muO_interp(x{1}, p);
    muG = @(p, T, y) 1./muG_interp(y{2}, p);
    rhoO = @(p, T, x) rhoO_interp(x{1}, p);
    rhoG = @(p, T, y) rhoG_interp(y{2}, p);
    
    isStableFn = [];
    mw = [1, 1];
    % mw = [model.fluid.rhoOS, model.fluid.rhoGS];

    cfluid = CompositionalFluid({'PseudoOil', 'PseudoGas'}, [nan, nan], [nan, nan], [nan, nan], [nan, nan], mw);

    eos = EquilibriumConstantModel(bomodel.G, cfluid, {K_o, K_g});
    prop = BlackOilPropertyModel({rhoO, rhoG}, {muO, muG}, isStableFn, eos.fluid);
    eos.PropertyModel = prop;
    
    f = bomodel.fluid;
    flds = {'rsSat', 'bO', 'bG', 'muO', 'muG'};
    for i = 1:numel(flds)
        if isfield(f, flds{i})
            f = rmfield(f, flds{i});
        end
    end
    model = NaturalVariablesCompositionalModel(bomodel.G, bomodel.rock, f, eos);

end