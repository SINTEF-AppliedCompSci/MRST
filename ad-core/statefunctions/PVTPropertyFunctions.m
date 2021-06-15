classdef PVTPropertyFunctions < StateFunctionGrouping
    % Default grouping for describing a system of flow equations. Contains
    % basic properties like mobility, density, component total mass etc.
    properties
        SurfaceDensity % Density at surface conditions (1 x nphase)
        Density % Phase density (1 x nphase)
        ShrinkageFactors % b = 1/B for each phase (1 x nphase)
        Viscosity % Phase viscosity (1 x nphase)
        PoreVolume % Effective pore-volumes (scalar)
        PhasePressures % Pressure of each phase (1 x nphase)
        PressureReductionFactors % Weighting factors to compute pressure equation (ncomp x 1)
    end

    methods
        function props = PVTPropertyFunctions(model)
            props@StateFunctionGrouping('PVTProps');
            pvt = props.getRegionPVT(model);

            % PVT properties
            isBO = isa(model, 'ThreePhaseBlackOilModel') && (model.disgas || model.vapoil);
            props = props.setStateFunction('Density', BlackOilDensity(model, pvt));
            if isBO
                mu = BlackOilViscosity(model, pvt);
                bfactors = BlackOilShrinkageFactors(model, pvt);
                % Black-oil specific features follow
                if model.disgas
                    props = props.setStateFunction('RsMax', RsMax(model, pvt));
                end
                % Vaporized oil
                if model.vapoil
                    props = props.setStateFunction('RvMax', RvMax(model, pvt));
                end
                w_p = BlackOilPressureReductionFactors(model);
            else
                mu = Viscosity(model, pvt);
                bfactors = ShrinkageFactors(model, pvt);
                w_p = PressureReductionFactors(model);
            end
            props = props.setStateFunction('SurfaceDensity', SurfaceDensity(model, pvt));
            props = props.setStateFunction('Viscosity', mu);
            props = props.setStateFunction('ShrinkageFactors', bfactors);
            props = props.setStateFunction('PressureReductionFactors', w_p);
            props = props.setStateFunction('PhasePressures', PhasePressures(model));
            if isfield(model.fluid, 'pvMultR')
                % Check for multiplier
                pv = BlackOilPoreVolume(model, pvt);
            else
                pv = PoreVolume(model, pvt);
            end
            props = props.setStateFunction('PoreVolume', pv);
            
        end
        
        function pvt = getRegionPVT(props, model)
            % Get the region indicator in each cell of the domain from the
            % rock. If not found in rock, region 1 is used in all cells.
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
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
