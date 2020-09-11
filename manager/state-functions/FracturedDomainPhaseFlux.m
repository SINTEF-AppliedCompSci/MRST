
classdef FracturedDomainPhaseFlux < StateFunction
    properties

    end
    
    methods
        function fm = FracturedDomainPhaseFlux(model)
            fm@StateFunction(model);
            fm = fm.dependsOn({'FaceMobility', 'PermeabilityPotentialGradient'});
        end

        
        function v = evaluateOnDomain(prop, model, state)
            
            [mob, kgrad] = prop.getEvaluatedDependencies(state,...
                'FaceMobility', 'PermeabilityPotentialGradient');
            
            %% Standard flux evaluation
            v = cellfun(@(x,y)-x.*y, mob, kgrad, 'UniformOutput', false);
            
            %% Compute fluxes between mutli-continuum domains
            for j = 1:length(model.G.FracturedDomains.domains)
                dom = model.G.FracturedDomains.domains{j};
                if(strcmp(dom.type,'multi_continuum'))
                    ids = dom.global_connection_ids;
                    vf = dom.transfer_model.transfer(model, state, j);
                    for i = 1:numel(mob)
                        v{i}(ids) = vf{i};
                    end
                end
            end
        end
    end
end