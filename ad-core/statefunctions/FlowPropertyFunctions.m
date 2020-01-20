classdef FlowPropertyFunctions < StateFunctionGrouping
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
        function props = FlowPropertyFunctions(model)
            props@StateFunctionGrouping();
            sat = props.getRegionSaturation(model);
            pvt = props.getRegionPVT(model);
            % Saturation properties
            props = props.setStateFunction('CapillaryPressure', BlackOilCapillaryPressure(model, sat));
            props = props.setStateFunction('RelativePermeability', BaseRelativePermeability(model, sat));
            props = props.setStateFunction('Mobility', Mobility(model, sat));

            % PVT properties
            props = props.setStateFunction('ShrinkageFactors', BlackOilShrinkageFactors(model, pvt));
            props = props.setStateFunction('Density', BlackOilDensity(model, pvt));
            props = props.setStateFunction('Viscosity', BlackOilViscosity(model, pvt));
            if isfield(model.fluid, 'pvMultR')
                % Check for multiplier
                pv = BlackOilPoreVolume(model, pvt);
            else
                pv = PoreVolume(model, pvt);
            end
            props = props.setStateFunction('PoreVolume', pv);
            props = props.setStateFunction('PhasePressures', PhasePressures(model, pvt));
            props = props.setStateFunction('PressureReductionFactors', BlackOilPressureReductionFactors(model));
            
            % Components
            props = props.setStateFunction('ComponentPhaseMass', ComponentPhaseMass(model));
            props = props.setStateFunction('ComponentTotalMass', ComponentTotalMass(model));
            props = props.setStateFunction('ComponentMobility', ComponentMobility(model));
            props = props.setStateFunction('ComponentPhaseDensity', ComponentPhaseDensity(model));

            if ~isempty(model.inputdata)
                % We may have recieved a deck. Check for endpoint scaling
                deck = model.inputdata;

                % Scaling is active and can impact rel.perm and pc
                do_scaling = isfield(deck.RUNSPEC, 'ENDSCALE');
                props.CapillaryPressure.scalingActive = do_scaling;
                props.RelativePermeability.scalingActive = do_scaling;

                % Set number of points for scaling
                three_point = isfield(deck.PROPS, 'SCALECRS') && strcmpi(deck.PROPS.SCALECRS{1}(1), 'y');
                props.RelativePermeability.relpermPoints = 2 + three_point;
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

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
