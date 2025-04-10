classdef FoamFlowPropertyFunctions < FlowPropertyFunctions
    
    properties
        ConcentrationsPartitioning
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function props = FoamFlowPropertyFunctions(model)
        % Retrieve necessary properties for the foam effect calculations
            
            % Parent constructor
            props = props@FlowPropertyFunctions(model);
            % Saturation regions
            satreg  = props.getRegionSaturation(model);
            surfreg = props.getRegionSurfactant(model); %#ok
            % Concentration partitioning into gas/water
            props = props.setStateFunction('ConcentrationsPartitioning', ...
                ConcentrationsPartitioning(model, satreg));
            % Mobility with foam-dependent multiplier for gas phase
            props = props.setStateFunction('Mobility', ...
                MultipliedMobility(model, satreg));
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function sat = getRegionSurfactant(props, model)
            r = model.rock;
            sat = ones(model.G.cells.num, 1);
            if isfield(r, 'regions')
                if isfield(r.regions, 'surfactant')
                    sat = r.regions.surfactant;
                end
            end
        end
        %-----------------------------------------------------------------%
        
    end
    
end