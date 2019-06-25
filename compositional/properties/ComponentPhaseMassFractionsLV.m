classdef ComponentPhaseMassFractionsLV < StateFunction
    properties
    end
    
    methods
        function gp = ComponentPhaseMassFractionsLV(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'ComponentPhaseMoleFractions'});
        end

        function mass = evaluateOnDomain(prop, model, state)
            eos = model.EOSModel;
            moles = prop.getEvaluatedDependencies(state, 'ComponentPhaseMoleFractions');
            
            mass = moles;
            offset = model.water + 1;
            for i = 1:size(mass, 2)
                mass(offset:end, i) = eos.getMassFraction(moles(offset:end, i));
            end
        end
    end
end