classdef PressureGradient < GridProperty
    properties

    end
    
    methods
        function gp = PressureGradient(varargin)
            gp@GridProperty(varargin{:});
        end
        function dp = evaluateOnDomain(prop, model, state)
            act = model.getActivePhases();
            nph = sum(act);
            Grad = model.operators.Grad;
            
            dp = cell(1, nph);
            if model.FlowPropertyFunctions.CapillaryPressure.pcPresent(model)
                % We have different phase pressures, call gradient once for
                % each phase
                p = model.getProp(state, 'PhasePressures');
                for i = 1:nph
                    dp{i} = Grad(p{i});
                end
            else
                % There is no capillary pressure and a single gradient for
                % the unique pressure is sufficient
                p = model.getProp(state, 'pressure');
                [dp{:}] = deal(Grad(p));
            end
        end
    end
end