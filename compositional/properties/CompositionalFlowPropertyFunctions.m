classdef CompositionalFlowPropertyFunctions < FlowPropertyFunctions
    properties
        ComponentPhaseMassFractions
        ComponentPhaseMoleFractions
    end
    
    methods
        function props = CompositionalFlowPropertyFunctions(model)
            props@FlowPropertyFunctions(model);
            pvt = props.getRegionPVT(model);
            props = props.setStateFunction('ShrinkageFactors', DensityDerivedShrinkageFactors(model, pvt));
            props = props.setStateFunction('Density', CompositionalDensity(model, pvt));
            
            props = props.setStateFunction('ComponentPhaseMassFractions', ComponentPhaseMassFractionsLV(model));
            props = props.setStateFunction('ComponentPhaseMoleFractions', ComponentPhaseMoleFractionsLV(model));

            props = props.setStateFunction('PhaseMixingCoefficients', PhaseMixingCoefficientsLV(model));
            props = props.setStateFunction('Fugacity', FugacityLV(model));
            props = props.setStateFunction('PhaseCompressibilityFactors', PhaseCompressibilityFactorsLV(model));
            props = props.setStateFunction('Viscosity', CompositionalViscosityLV(model));
        end
        
        function props = setCompactEvaluation(props, val)
            names = props.functionNames;
            for i = 1:numel(names)
                name = names{i};
                fn = props.getStateFunction(name);
                if isprop(fn, 'useCompactEvaluation')
                    fn.useCompactEvaluation = val;
                    props = props.setStateFunction(name, fn);
                end
            end
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
