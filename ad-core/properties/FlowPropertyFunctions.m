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
            % Saturation properties
            props.CapillaryPressure = BlackOilCapillaryPressure(model.AutoDiffBackend, sat);
            props.RelativePermeability = BaseRelativePermeability(model.AutoDiffBackend, sat);
            props.Mobility = Mobility(model.AutoDiffBackend, sat);

            % PVT properties
            props.ShrinkageFactors = BlackOilShrinkageFactors(model.AutoDiffBackend, pvt);
            props.Density = BlackOilDensity(model.AutoDiffBackend, pvt);
            props.Viscosity = BlackOilViscosity(model.AutoDiffBackend, pvt);
            props.PoreVolume = MultipliedPoreVolume(model.AutoDiffBackend, pvt);
            props.PhasePressures = PhasePressures(model.AutoDiffBackend, pvt);
            if ~isempty(model.inputdata)
                deck = model.inputdata;
                
                do_scaling = isfield(deck.RUNSPEC, 'ENDSCALE');
                three_point = isfield(deck.PROPS, 'SCALECRS') && strcmpi(deck.PROPS.SCALECRS{1}(1), 'y');
                props.RelativePermeability.relpermScaling = do_scaling;
                props.RelativePermeability.relpermPoints = 2 + three_point;
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