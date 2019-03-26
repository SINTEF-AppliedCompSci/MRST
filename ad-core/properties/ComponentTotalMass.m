classdef ComponentTotalMass <  GridProperty & ComponentProperty
    properties

    end
    
    methods
        function gp = ComponentTotalMass(model, varargin)
            gp@GridProperty(model, varargin{:});
            gp = gp.dependsOn('ComponentPhaseMass');
        end
        function mass = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            mass = cell(ncomp, 1);
            phase_mass = state.(prop.structName).ComponentPhaseMass;
            for c = 1:ncomp
                % Loop over phases where the component may be present
                for ph = 1:nph
                    % Check if present
                    m = phase_mass{c, ph};
                    if ~isempty(m)
                        if isempty(mass{c})
                            mass{c} = m;
                        else
                            mass{c} = mass{c} + m;
                        end
                    end
                end
            end
        end
    end
end