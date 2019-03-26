classdef FlowPropertyFunctions < PropertyFunctions
    properties
        Density
        Viscosity
        RelativePermeability
        CapillaryPressure
        PhasePressures
        ShrinkageFactors
        Mobility
        PoreVolume
        
        ComponentTotalMass % Total component mass
        ComponentPhaseMass % Component mass in each phase
        ComponentMobility
        ComponentPhaseDensity
    end
    
    methods
        function props = FlowPropertyFunctions(model)
            r = model.rock;
            props@PropertyFunctions();
            [sat, pvt] = deal(ones(model.G.cells.num, 1));
            if isfield(r, 'regions')
                if isfield(r.regions, 'saturation')
                    sat = r.regions.saturation;
                end
                if isfield(r.regions, 'pvt')
                    pvt = r.regions.pvt;
                end
            end
            ad = model.AutoDiffBackend;
            % Saturation properties
            props.CapillaryPressure = BlackOilCapillaryPressure(model, sat);
            props.RelativePermeability = BaseRelativePermeability(model, sat);
            props.Mobility = Mobility(model, sat);

            % PVT properties
            props.ShrinkageFactors = BlackOilShrinkageFactors(model, pvt);
            props.Density = BlackOilDensity(model, pvt);
            props.Viscosity = BlackOilViscosity(model, pvt);
            props.PoreVolume = MultipliedPoreVolume(model, pvt);
            props.PhasePressures = PhasePressures(model, pvt);
            
            % Components
            props.ComponentPhaseMass = ComponentPhaseMass(model);
            props.ComponentTotalMass = ComponentTotalMass(model);
            props.ComponentMobility = ComponentMobility(model);
            props.ComponentPhaseDensity = ComponentPhaseDensity(model);

            if ~isempty(model.inputdata)
                deck = model.inputdata;
                
                do_scaling = isfield(deck.RUNSPEC, 'ENDSCALE');
                three_point = isfield(deck.PROPS, 'SCALECRS') && strcmpi(deck.PROPS.SCALECRS{1}(1), 'y');
                props.RelativePermeability.relpermScaling = do_scaling;
                props.RelativePermeability.relpermPoints = 2 + three_point;
                if isfield(model.rock, 'sw')
                    % Endpoint capillary pressure is defined
                    props.CapillaryPressure =...
                        props.CapillaryPressure.setWaterEndpointScaling...
                                                (model, model.rock.sw, 1);
                end
            end
            
            % Define storage
            props.structName = 'FlowProps';
        end

        function state = evaluateProperty(props, model, state, name)
            switch name
                case 'Mobility'
                    state = props.evaluateDependencies(model, state, {'Viscosity', 'RelativePermeability'});
                case {'Viscosity', 'ShrinkageFactors'}
                    state = props.evaluateDependencies(model, state, {'CapillaryPressure'});
                case {'Density'}
                    state = props.evaluateDependencies(model, state, {'ShrinkageFactors'});
            end
            state = evaluateProperty@PropertyFunctions(props, model, state, name);
        end
    end
end