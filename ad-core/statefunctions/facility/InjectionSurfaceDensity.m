classdef InjectionSurfaceDensity < StateFunction
    % Get injection surface density
    properties

    end
    
    methods
        function gp = InjectionSurfaceDensity(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('FacilityWellMapping');
            gp = gp.dependsOn('SurfaceDensity', 'PVTPropertyFunctions');
            gp.label = '\rho_\alpha^{w}';
        end
        function rhoS = evaluateOnDomain(prop, facility, state)
            model = facility.ReservoirModel;
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            % rhoSr = prop.getEvaluatedExternals(model, state, 'SurfaceDensity');
            rhoSr = model.PVTPropertyFunctions.get(model, state, 'SurfaceDensity');
            W = map.W;
            % We take the surface density for the first well cell,
            % regardless of active or inactive status for that
            % perforation.
            topcell = arrayfun(@(x) x.cells(1), W);
            rhoS = applyFunction(@(x) x(topcell), rhoSr);
            isInj = map.isInjector;
            if isfield(W, 'rhoS') && any(isInj)
                % Surface density is given on a per-well-basis for the
                % injectors
                rhoS_inj = vertcat(W(isInj).rhoS);
                for i = 1:numel(rhoS)
                    rhoS{i}(isInj) = rhoS_inj(:, i);
                end
            end
        end
    end
end

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
