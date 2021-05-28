classdef PolymerPermReduction < StateFunction
% Polymer permeability reduction factor. This coefficient divides the water
% relative permeability and accounts for effect the reduction of permeability
% that adsorped polymer induces.
    
    properties
    end

    methods
        function gp = PolymerPermReduction(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn('PolymerAdsorption');
            assert(model.water);
            assert(all(isfield(model.fluid, {'rrf','adsMax'})));
        end

        function permRed = evaluateOnDomain(prop, model, state)
            ads = prop.getEvaluatedDependencies(state,'PolymerAdsorption');
            fluid = model.fluid;
            permRed = 1 + ((fluid.rrf - 1)./fluid.adsMax).*ads;
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
