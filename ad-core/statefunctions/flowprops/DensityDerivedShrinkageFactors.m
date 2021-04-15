classdef DensityDerivedShrinkageFactors < StateFunction
    % Simple "shrinkage factors" which are just phase densities divided by
    % surface densities. Useful for models that provide density directly.
    properties
    end
    
    methods
        function gp = DensityDerivedShrinkageFactors(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'Density'});
            gp.label = 'b_\alpha';
        end

        function b = evaluateOnDomain(prop, model, state)
            b = prop.getEvaluatedDependencies(state, 'Density');
            rhoS = model.getSurfaceDensities();
            for i = 1:numel(b)
                b{i} = b{i}./rhoS(i);
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
