classdef ComponentImplementation
    properties
        name
    end
    
    methods
        function c = ComponentImplementation(name)
            c.name = name;
        end
        
        function [c, phasenames] = getComponentDensity(component, model, state, phasenames)
            if nargin < 4
                phasenames = model.getPhaseNames();
            end
            nph = numel(phasenames);
            c = cell(nph, 1);
        end

        function c = getPhaseComposition(component, model, state, varargin)
            c = component.getComponentDensity(model, state, varargin{:});
            rho = model.getProp(state, 'Density');
            c = cellfun(@rdivide, c, rho, 'UniformOutput', false);
        end
        
        function mass = getComponentMass(component, model, state, varargin)
            c = component.getComponentDensity(model, state, varargin{:});
            pv = model.getProp(state, 'PoreVolume');
            mass = cellfun(@(x) x.*pv, c, 'UniformOutput', false);
        end
        
        function mass = getComponentMobility(component, model, state, varargin)
            c = component.getComponentDensity(model, state, varargin{:});
            mob = model.getProp(state, 'Mobility');
            mass = cellfun(@(x) x.*mob, c, 'UniformOutput', false);
        end
    end
end