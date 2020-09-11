classdef FracturedDomainFlowPropertyFunctions < StateFunctionGrouping
    % Default grouping for describing a system of flow equations. Contains
    % basic properties like mobility, density, component total mass etc.
    properties
        Density % Phase density (1 x nphase)
        ShrinkageFactors % b = 1/B for each phase (1 x nphase)
        Viscosity % Phase viscosity (1 x nphase)
        RelativePermeability % Phase relative permeabilities (1 x nphase)
        CapillaryPressure % Capillary pressures (if present) (1 x nphase, with empty entries for phases with equal pressure to reference)
        PhasePressures % Phase pressures for each phase (1 x nphase)
        Mobility % Phase mobilities (1 x nphase)
        PoreVolume % Effective pore-volumes (scalar)
        PressureReductionFactors % Weighting factors to compute pressure equation (ncomp x 1)
        % Component properties
        ComponentTotalMass % Total component mass (ncomp x 1)
        ComponentPhaseMass % Component mass in each phase (ncomp x nphase)
        ComponentMobility % Mobility of each component in each phase (mass, not volume) (ncomp x nphase)
        ComponentPhaseDensity % Mass-density per volume in each phase (ncomp x nphase)
    end

    methods
        function props = FracturedDomainFlowPropertyFunctions(model)
            props@StateFunctionGrouping();
            sat = props.getRegionSaturation(model);
            pvt = props.getRegionPVT(model);
            % Saturation properties
            props.CapillaryPressure = BlackOilCapillaryPressure(model, sat);
            props.RelativePermeability = BaseRelativePermeability(model, sat);
            props.Mobility = Mobility(model, sat);

            % PVT properties
            props.ShrinkageFactors = BlackOilShrinkageFactors(model, pvt);
            props.Density = BlackOilDensity(model, pvt);
            props.Viscosity = BlackOilViscosity(model, pvt);
            props.PoreVolume = FracturedDomainPoreVolume(model, pvt);
            props.PhasePressures = PhasePressures(model, pvt);
            props.PressureReductionFactors = BlackOilPressureReductionFactors(model);
            
            % Components
            props.ComponentPhaseMass = ComponentPhaseMass(model);
            props.ComponentTotalMass = ComponentTotalMass(model);
            props.ComponentMobility = ComponentMobility(model);
            props.ComponentPhaseDensity = ComponentPhaseDensity(model);

            if ~isempty(model.inputdata)
                % We may have recieved a deck. Check for endpoint scaling
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
            % Black-oil specific features follow
            if isprop(model, 'disgas') && model.disgas
                props = props.setStateFunction('RsMax', RsMax(model, pvt));
            end
            % Vaporized oil
            if isprop(model, 'vapoil') && model.vapoil
                props = props.setStateFunction('RvMax', RvMax(model, pvt));
            end
            % Define storage field in state
            props.structName = 'FlowProps';
        end
        
        function sat = getRegionSaturation(props, model)
            r = model.rock;
            sat = ones(model.G.cells.num, 1);
            if isfield(r, 'regions')
                if isfield(r.regions, 'saturation')
                    sat = r.regions.saturation;
                end
            end
        end
        
        function pvt = getRegionPVT(props, model)
            r = model.rock;
            pvt = ones(model.G.cells.num, 1);
            if isfield(r, 'regions')
                if isfield(r.regions, 'pvt')
                    pvt = r.regions.pvt;
                end
            end
        end
    end
end