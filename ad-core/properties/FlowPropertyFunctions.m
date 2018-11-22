classdef FlowPropertyFunctions < PropertyFunctions
    properties
        Density
        Viscosity
        RelativePermeability
        CapillaryPressure
    end
    
    methods
        function props = FlowPropertyFunctions(model)
            r = model.rock;
            [sat, pvt] = deal([]);
            if isfield(r, 'regions')
                if isfield(r.regions, 'saturation')
                    sat = r.regions.saturation;
                end
                if isfield(r.regions, 'pvt')
                    pvt = r.regions.pvt;
                end
            end
            props.Density = BlackOilDensity(model.AutoDiffBackend, pvt);
        end
    end
end