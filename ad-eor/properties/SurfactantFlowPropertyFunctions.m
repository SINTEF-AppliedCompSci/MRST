classdef SurfactantFlowPropertyFunctions < FlowPropertyFunctions
    properties
        CapillaryNumber
        SurfactantAdsorption
    end
    
    methods
        function props = SurfactantFlowPropertyFunctions(model)
            props = props@FlowPropertyFunctions(model);
            satreg  = props.getRegionSaturation(model);
            surfreg = props.getRegionSurfactant(model);
            
            props.RelativePermeability = SurfactantRelativePermeability(model, ...
                                                              satreg, surfreg);
            props.CapillaryPressure    = SurfactantCapillaryPressure(model, satreg);
            props.CapillaryNumber      = CapillaryNumber(model);
            props.SurfactantAdsorption = SurfactantAdsorption(model);
            props.Viscosity = BlackOilSurfactantViscosity(model);
        end
        
        function sat = getRegionSurfactant(props, model)
            r = model.rock;
            sat = ones(model.G.cells.num, 1);
            if isfield(r, 'regions')
                if isfield(r.regions, 'surfactant')
                    sat = r.regions.surfactant;
                end
            end
        end
        
    end
end