classdef ComponentImplementation
    % Base class for all component class instances
    properties
        name
        dependencies = {};
        externals = [];
    end
    
    methods
        function c = ComponentImplementation(name)
            c.name = name;
            % Document dependencies internal to grouping
            c = c.dependsOn({'PoreVolume', 'Density', 'Mobility'});
            % State dependencies
            c = c.dependsOn('s', 'state');
        end
        
        function c = getComponentDensity(component, model, state)
            % Density of component in each phase
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
        end

        function c = getPhaseCompositionWell(component, model, state, W)
            % Get composition of phases in well (on injection)
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
        end
        
        function c = getPhaseCompositionSurface(component, model, state, pressure, temperature)
            % Surface compositon
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
        end
        
        function c = getPurePhaseDensitySurface(component, model, state, pressure, temperature)
            % Surface density, for a pure component
            nph = model.getNumberOfPhases();
            c = cell(nph, 1);
        end

        function c = getPhaseComposition(component, model, state)
            % Density of component in each phase. Default implementation
            % goes via the component density and total density and assumes
            % these are internally consistent.
            c = component.getComponentDensity(model);
            rho = model.getProp(state, 'Density');
            c = cellfun(@rdivide, c, rho, 'UniformOutput', false);
        end
        
        function mass = getComponentMass(component, model, state, varargin)
            % Mass of component in each phase
            % saturation * pore-volume * component mass density
            pv = model.getProp(state, 'PoreVolume');
            mass = component.getComponentDensity(model, state, varargin{:});
            ph = model.getPhaseNames();
            % Iterate over phases and weight by pore-volume and saturation
            for i = 1:numel(mass)
                if ~isempty(mass{i})
                    s = model.getProp(state, ['s', ph(i)]);
                    mass{i} = s.*pv.*mass{i};
                end
            end
        end
        
        function cmob = getComponentMobility(component, model, state, varargin)
            % Product of mobility and component phase density
            mass = component.getComponentDensity(model, state, varargin{:});
            mob = model.getProp(state, 'Mobility');
            
            nphase = numel(mass);
            cmob = cell(1, nphase);
            for i = 1:nphase
                if ~isempty(mass{i})
                    cmob{i} = mob{i}.*mass{i};
                end
            end
        end

        function prop = dependsOn(prop, varargin)
            % Document dependencies and external dependencies
            prop = addPropertyDependence(prop, varargin{:});
        end
    end
end