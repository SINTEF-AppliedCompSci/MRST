function [eqs, names, types, state] = biotCompositionalEquations(model, state0, state, dt, drivingForces, varargin)
%Undocumented Utility Function

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


    G     = model.G;
    fluid = model.fluid;
    op    = model.operators;

    pv         = op.pv;
    momentop   = op.momentop;
    mechdirop  = op.mechDirichletop;

    u  = model.getProp(state, 'u');
    bp = model.getProp(state, 'bp');
    lm = model.getProp(state, 'lambdamech');

    u0  = model.getProp(state0, 'u');
    bp0 = model.getProp(state0, 'bp');
    lm0 = model.getProp(state0, 'lambdamech');

    p = model.getProp(state, 'pressure');

    fac = 1 / (1e6 * mean(G.cells.volumes));
    meqs{1} = fac*momentop(u, p, lm);
    meqs{2} = fac*mechdirop(u, p, lm);
    meqs{3} = p - bp;

    mnames = {'momentum', 'mechbcs', 'coupling'};
    mtypes  = {'cells', 'bc', 'cells'};

    [feqs, flux, fnames, ftypes] = model.FlowDiscretization.componentConservationEquations(model, state, state0, dt);
    src = model.FacilityModel.getComponentSources(state);

    % Treat source or bc terms
    if ~isempty(drivingForces.bc) || ~isempty(drivingForces.src)
        [pressures, sat, mob, rho, rs, rv] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'Rs', 'Rv');
        dissolved = model.getDissolutionMatrix(rs, rv);
        feqs = model.addBoundaryConditionsAndSources(feqs, fnames, ftypes, state, ...
                                                    pressures, sat, mob, rho, ...
                                                    dissolved, {}, ...
                                                    drivingForces);
    end

    % Add aquifer contributions if any.
    if ~isempty(model.AquiferModel)
        feqs = addAquifersContribution(model.AquiferModel, feqs, fnames, state, dt);
    end

    % Add sources
    feqs = model.insertSources(feqs, src);
    % Assemble equations
    for i = 1:numel(feqs)
        feqs{i} = model.operators.AccDiv(feqs{i}, flux{i});
    end

    % Get facility equations
    [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);

    eqs = [meqs, feqs, weqs];
    names = [mnames, fnames, wnames];
    types = [mtypes, ftypes, wtypes];
end
