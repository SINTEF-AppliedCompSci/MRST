classdef ComponentPhaseFlux < StateFunction
    properties (Access = protected)
        mobility_name; % Name of state function where component mobility comes from
    end

    methods
        function cf = ComponentPhaseFlux(model, mob_name)
            cf@StateFunction(model);
            if nargin < 2
                mob_name = 'FaceComponentMobility'; 
            end
            cf.mobility_name = mob_name;
            cf = cf.dependsOn({'PermeabilityPotentialGradient', cf.mobility_name});
        end

        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            [kgrad, compMob] = prop.getEvaluatedDependencies(state,...
                'PermeabilityPotentialGradient', prop.mobility_name);
            v = cell(ncomp, nph);
            for c = 1:ncomp
                for ph = 1:nph
                    mob = compMob{c, ph};
                    if ~isempty(mob)
                        v{c, ph} = -mob.*kgrad{ph};
                    end
                end
            end
        end
    end
end