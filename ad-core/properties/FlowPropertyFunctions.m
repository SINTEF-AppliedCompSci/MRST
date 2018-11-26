classdef FlowPropertyFunctions < PropertyFunctions
    properties
        Density
        Viscosity
        RelativePermeability
        CapillaryPressure
        ShrinkageFactors
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
            props.ShrinkageFactors = BlackOilShrinkageFactors(model.AutoDiffBackend, pvt);
            props.Density = BlackOilDensity(model.AutoDiffBackend, pvt);
            props.Viscosity = BlackOilViscosity(model.AutoDiffBackend, pvt);
            
            if ~isempty(model.inputdata)
                deck = model.inputdata;
                
                do_scaling = isfield(deck.RUNSPEC, 'ENDSCALE');
                three_point = isfield(deck.PROPS, 'SCALECRS') && strcmpi(deck.PROPS.SCALECRS{1}(1), 'y');
                props.RelativePermeability.relpermScaling = do_scaling;
                props.RelativePermeability.relpermPoints = 2 + three_point;
            end
            
            % Define storage
            props.structName = 'FlowProps';
            props.structFields = {'CapillaryPressure', ...
                                  'Density', ...
                                  'RelativePermeability', ...
                                  'ShrinkageFactors', ...
                                  'Viscosity', ...
                                  'Mobility', ...
                                  'Components'};
        end

        function evaluateProperty(props, model, state, name)
            switch name
                case 'Mobility'
                    props.evaluateDependencies(model, state, {'Viscosity', 'RelativePermeability'});
                case {'Viscosity', 'ShrinkageFactors'}
                    props.evaluateDependencies(model, state, {'CapillaryPressure'});
                case {'Density'}
                    props.evaluateDependencies(model, state, {'ShrinkageFactors'});
            end
            evaluateProperty@PropertyFunctions(props, model, state, name);
        end
    end
end