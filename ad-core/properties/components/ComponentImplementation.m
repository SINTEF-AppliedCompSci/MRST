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
            pv = model.getProp(state, 'PoreVolume');
            mass = component.getComponentDensity(model, state, varargin{:});
            for i = 1:numel(mass)
                if ~isempty(mass{i})
                    mass{i} = pv.*mass{i};
                end
            end
        end
        
        function cmob = getComponentMobility(component, model, state, varargin)
            mass = component.getComponentDensity(model, state, varargin{:});
            mob = model.getProp(state, 'Mobility');
            
            ncomp = numel(mass);
            cmob = cell(1, ncomp);
            for i = 1:ncomp
                if ~isempty(mass{i})
                    cmob{i} = mob{i}.*mass{i};
                end
            end
        end

    end
end