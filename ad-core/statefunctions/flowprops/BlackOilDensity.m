classdef BlackOilDensity < StateFunction
    properties
        disgas = false;
        vapoil = false;
    end
    methods
        function gp = BlackOilDensity(model, varargin)
            gp@StateFunction(model, varargin{:});
            if isprop(model, 'disgas')
                gp.disgas = model.disgas;
                if gp.disgas
                    gp = gp.dependsOn({'rs'}, 'state');
                end
            end
            if isprop(model, 'vapoil')
                gp.vapoil = model.vapoil;
                if gp.vapoil
                    gp = gp.dependsOn({'rv'}, 'state');
                end
            end
            gp = gp.dependsOn({'ShrinkageFactors'});
            gp.label = '\rho_\alpha';
        end
        function rho = evaluateOnDomain(prop, model, state)
            rhoS = model.getSurfaceDensities(prop.regions);
            nph = size(rhoS, 2);
            rho = cell(1, nph);
            b = prop.getEvaluatedDependencies(state, 'ShrinkageFactors');
            for i = 1:nph
                rho{i} = rhoS(:, i).*b{i};
            end
            if (prop.disgas || prop.vapoil) && model.gas && model.oil
                [oix, gix] = model.getPhaseIndex('O', 'G');
                if prop.disgas
                    rs = model.getProp(state, 'rs');
                    rho{oix} = rho{oix} + rs.*b{oix}.*rhoS(:, gix);
                end
                if prop.vapoil
                    rv = model.getProp(state, 'rv');
                    rho{gix} = rho{gix} + rv.*b{gix}.*rhoS(:, oix);
                end
            end
            mv = cellfun(@(x) min(value(x)), rho);
            if any(mv <= 0)
                if model.verbose > 1
                    warning('Negative densities detected! Capping to 1e-12.')
                end
                rho = cellfun(@(x) max(x, 1e-12), rho, 'UniformOutput', false);
            end
        end
    end
end

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
