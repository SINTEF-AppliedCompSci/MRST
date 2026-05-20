function forces = validateCompositionalForcesSW(model, forces, it)
% Soreide-Whitson specific extension of validateCompositionalForces.
% Computes well compi using flash-based saturations for the SW EOS.
%
% This is called after the standard validateCompositionalForces and
% overrides the compi values with flash-based saturations, which is
% needed for water-containing mixtures with the Soreide-Whitson EOS.

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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
        T = model.FacilityModel.T;
        p = model.FacilityModel.pressure;
        eos = model.EOSModel;
        act = vertcat(forces.W.status);
        if any(act)
            wellIndices = find(act);
            z = vertcat(forces.W(act).components);
            n = size(z, 1);
            [L, x, y, Z_L, Z_V, rhoL, rhoV] = standaloneFlash(repmat(p, n, 1), repmat(T, n, 1), z, eos);
            for i = 1:numel(wellIndices)
                wNo = wellIndices(i);
                [sL, sV] = eos.computeSaturations(p, T, rhoL(i), rhoV(i), x(i, :), y(i, :), L(i), Z_L(i), Z_V(i));
                L_mass = sL.*rhoL(i)./(sL.*rhoL(i) + sV.*rhoV(i));
                comp = [L_mass, 1-L_mass];
                forces.W(wNo).compi = comp;
                forces.W(wNo).rhoS = model.getSurfaceDensities();
                if isfield(model.rock, 'regions')
                    if isfield(model.rock.regions, 'pvt')
                        forces.W(wNo).rhoS = forces.W(wNo).rhoS(model.rock.regions.pvt(forces.W(wNo).cells(1)), :);
                    end
                end
                forces.W(wNo).rhoS(model.getLiquidIndex()) = rhoL(i);
                forces.W(wNo).rhoS(model.getVaporIndex()) = rhoV(i);
            end
        end
    end
end
