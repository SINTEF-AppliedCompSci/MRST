function model = convertBlackOilModelToCompositionalModel(bomodel, varargin)
% Convert blackoil model to compositional-like model

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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

    cfluid = CompositionalMixture({'PseudoOil', 'PseudoGas'}, [nan, nan], [nan, nan], [nan, nan], [nan, nan], mw);

    eos = EquilibriumConstantModel(bomodel.G, cfluid, {K_o, K_g});
    prop = BlackOilPropertyModel({rhoO, rhoG}, {muO, muG}, isStableFn, eos.CompositionalMixture);
    eos.PropertyModel = prop;
    
    f = bomodel.fluid;
    flds = {'rsSat', 'bO', 'bG', 'muO', 'muG'};
    for i = 1:numel(flds)
        if isfield(f, flds{i})
            f = rmfield(f, flds{i});
        end
    end
    model = NaturalVariablesCompositionalModel(bomodel.G, bomodel.rock, f, eos);
    model.operators = bomodel.operators;
end
