classdef ImmiscibleComponent < ComponentImplementation
    properties
        phaseIndex % Index of phase this component belongs to
    end
    
    methods
        function c = ImmiscibleComponent(name, phase)
            c@ComponentImplementation(name);
            c.phaseIndex = phase;
        end
        
        function [c, phasenames] = getComponentDensity(component, model, state, varargin)
            [c, phasenames] = getComponentDensity@ComponentImplementation(component, model, state, varargin{:});
            c{component.phaseIndex} = 1;
        end
    end
end