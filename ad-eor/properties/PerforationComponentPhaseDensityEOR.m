classdef PerforationComponentPhaseDensityEOR < PerforationComponentPhaseDensity
    % Component density to used for each well connection
    properties

    end
    
    methods
        function gp = PerforationComponentPhaseDensityEOR(varargin)
            gp = gp@PerforationComponentPhaseDensity(varargin{:});
            % TODO: right the follow?
            gp = gp.dependsOn({'PolymerEffViscMult', 'PolymerViscMult'}, 'PVTPropertyFunctions');
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp.label = '\rho_{wc}';
        end
        
        
        function rhoc = evaluateOnDomain(prop, model, state)
            rhoc = evaluateOnDomain@PerforationComponentPhaseDensity(prop, model, state);
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            for c = 1 : model.getNumberOfComponents()
                comp = model.ReservoirModel.Components{c};
                if strcmp(comp.name, 'polymer')
                    wIx = 1;
                    % model.ReservoirModel
                    [effviscmult, pviscmult] = model.ReservoirModel.getProps(state, 'PolymerEffViscMult', 'PolymerViscMult');
                    effvismultw = effviscmult(map.cells);
                    pviscmultw = pviscmult(map.cells);
                    rhoc{c, wIx} = rhoc{c, wIx} .* effvismultw ./ pviscmultw;
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
