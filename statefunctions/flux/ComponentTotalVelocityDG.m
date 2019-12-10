classdef ComponentTotalVelocityDG < StateFunction
    properties
    end
    
    methods
        function gp = ComponentTotalVelocityDG(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('ComponentPhaseVelocity');
        end
        
        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();
            v = cell(ncomp, 1);
            phase_flux = prop.getEvaluatedDependencies(state, 'ComponentPhaseVelocity');
            
            for c = 1:ncomp
                % Loop over phases where the component may be present
                for ph = 1:nph
                    % Check if present
                    m = phase_flux{c, ph};
                    if ~isempty(m)
                        if isempty(v{c})
                            v{c} = m;
                        else
                            v{c} = v{c} + m;
                        end
                    end
                end
            end
        end
    end
end