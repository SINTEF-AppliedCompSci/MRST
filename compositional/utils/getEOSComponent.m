function c = getEOSComponent(model, p, T, name, ci)
%Undocumented Utility Function

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
    mixture = model.EOSModel.CompositionalMixture;
    names_hc = mixture.names;
    hcpos = strcmp(names_hc, name);
    n_hc = numel(names_hc);
    z = zeros(1, n_hc);
    z(hcpos) = 1;
    [L, ~, ~, ~, ~, rhoL, rhoV] = standaloneFlash(p, T, z, model.EOSModel);
    Lm = L.*rhoL./(rhoL.*L + rhoV.*(1-L));
    Vm = 1 - Lm;
    % Get indices
    nph = model.getNumberOfPhases();
    Li = model.getLiquidIndex();
    Vi = model.getVaporIndex();
    % Assign values
    frac = zeros(1, nph);
    frac(Li) = Lm;
    rho(Li) = rhoL;
    
    frac(Vi) = Vm;
    rho(Vi) = rhoV;
    extra = model.getNonEoSPhaseNames();
    phases = model.getPhaseNames();
    for i = 1:numel(extra)
        rho(phases == extra(i)) = model.fluid.(['rho', extra(i), 'S']);
    end
    c = EquationOfStateComponent(name, p, T, ci, frac, rho, mixture.molarMass(hcpos));
end
