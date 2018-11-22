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
            f = model.fluid;
            [sat, pvt] = deal(ones(model.G.cells.num, 1));
            if isfield(r, 'regions')
                if isfield(r.regions, 'saturation')
                    sat = r.regions.saturation;
                end
                if isfield(r.regions, 'pvt')
                    pvt = r.regions.pvt;
                end
            end
            % Saturation properties
            props.CapillaryPressure = BlackOilCapillaryPressure(model.AutoDiffBackend, sat);
            props.RelativePermeability = StoneRelativePermeability1(model.AutoDiffBackend, sat);
            
            % PVT properties
            props.Density = BlackOilDensity(model.AutoDiffBackend, pvt);
            props.Viscosity = BlackOilViscosity(model.AutoDiffBackend, pvt);
            
            % Define storage
            props.structName = 'FlowProps';
            props.structFields = {'CapillaryPressure', ...
                                  'Density', ...
                                  'RelativePermeability', ...
                                  'Viscosity', ...
                                  'Mobility', ...
                                  'Components'};
        end

        function evaluateProperty(props, model, state, name)
            switch name
                case 'Mobility'
                    props.evaluateDependencies(model, state, {'Viscosity', 'RelativePermeability'});
                case {'Density', 'Viscosity'}
                    props.evaluateDependencies(model, state, {'CapillaryPressure'});
            end
            evaluateProperty@PropertyFunctions(props, model, state, name);
        end
    end
end