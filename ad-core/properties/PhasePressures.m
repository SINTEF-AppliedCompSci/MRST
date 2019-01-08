classdef PhasePressures < GridProperty
    properties
    end
    
    methods
        function p_phase = evaluateOnGrid(prop, model, state)
            [pc, p] = model.getProps(state, 'CapillaryPressure', 'Pressure');
            
            nph = numel(pc);
            p_phase = cell(1, nph);
            for i = 1:nph
                if isempty(pc{i})
                    p_phase{i} = p;
                else
                    p_phase{i} = p + pc{i};
                end
            end
        end
    end
end