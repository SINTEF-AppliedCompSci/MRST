function forces = validateCompositionalForces(model, forces, it)
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

    if ~isempty(forces.W) && numel(forces.W) > 0
        assert(~isempty(model.FacilityModel), ...
            'FacilityModel must be set up before validating driving forces for a compositional problem with wells!');
        T = model.FacilityModel.T;
        p = model.FacilityModel.pressure;
        
        assert(isfield(forces.W, 'components'), ...
            'Wells must have field .components for a compositional model.');
        eos = model.EOSModel;
        if ~isfield(forces.W, 'rhoS')
            act = vertcat(forces.W.status);
            if any(act)
                wellIndices = find(act);
                z = vertcat(forces.W(act).components);
                Z = eos.getMassFraction(z);
                n = size(z, 1);
                [L, x, y, Z_L, Z_V, rhoL, rhoV] = standaloneFlash(repmat(p, n, 1), repmat(T, n, 1), z, eos);

                for i = 1:numel(wellIndices)
                    wNo = wellIndices(i);
                    [rho, comp] = getSurfaceParameters(model, forces.W(wNo), rhoL(i), rhoV(i), x(i, :), y(i, :), L(i), Z_L(i), Z_V(i), Z(i, :));
                    forces.W(wNo).compi = comp;
                    forces.W(wNo).rhoS = rho;
                end
            end
        end
    end
end

function [rho, compi] = getSurfaceParameters(model, W, rhoL, rhoV, x, y, L, Z_L, Z_V, Z)
    hc = model.getEoSPhaseIndices();
    nph = model.getNumberOfPhases();
    if false
        % Use flash
        [sL, sV] = eos.computeSaturations(rhoL, rhoV, x, y, L, Z_L, Z_V);
        % compi is a mass-fraction in practice
        L_mass = sL.*rhoL(u)./(sL.*rhoL + sV.*rhoV);
        comp = [L_mass, 1-L_mass];
    else
        % Use the pre-computed definition of light/heavy
        % components to determine "compi"
        isEOS = cellfun(@(x) isa(x, 'EquationOfStateComponent'), model.Components);
        val = cellfun(@(x) x.surfacePhaseMassFractions, model.Components(isEOS), 'UniformOutput', false)';
        val = vertcat(val{:});
        val = val(:, hc);
        comp = sum(bsxfun(@times, val, Z'), 1);
    end
    rho = model.getSurfaceDensities();
    rho(model.getLiquidIndex()) = rhoL;
    rho(model.getVaporIndex()) = rhoV;
    if numel(hc) < nph
        assert(~isempty(W.compi), ...
            'W.compi must be present for compositional flow with immiscible phases.');
        tmp = sum(W.compi(hc));
        compi = W.compi;
        compi(hc) = tmp.*comp;
    else
        compi = comp;
    end
end
