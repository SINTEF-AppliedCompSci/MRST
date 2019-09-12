classdef PolymerPhaseFlux < StateFunction
    properties

    end
    
    methods
        function gp = PolymerPhaseFlux(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PhaseFlux', 'FaceConcentration'});
        end
        
        function v = evaluateOnDomain(prop, model, state)
            [phaseFlux, cf] = prop.getEvaluatedDependencies(state,...
                'PhaseFlux', 'FaceConcentration');
            nph = numel(phaseFlux);
            v = cell(1, nph+1);
            for i = 1:nph
                v{i} = deal(phaseFlux{i});
            end
            
            fluid = model.fluid;
            mixpar = fluid.mixPar;
            cbar   = cf/fluid.cmax;
            a = fluid.muWMult(fluid.cmax).^(1-mixpar);
            v{nph+1} = v{1}./(1+(1-a)*cbar).*cf;
        end
    end
end