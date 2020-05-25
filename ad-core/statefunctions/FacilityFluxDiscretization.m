classdef FacilityFluxDiscretization < StateFunctionGrouping
    properties
        PhaseFlux
        ComponentTotalFlux
        ComponentPhaseFlux
        PerforationPressureGradient
        WellIndex
        FacilityWellMapping
        InjectionSurfaceDensity
        ComponentPhaseFractionInjectors
        PerforationMobility
        PerforationComponentPhaseDensity
    end
    
    methods
        function group = FacilityFluxDiscretization(model)
            group@StateFunctionGrouping('FacilityFluxProps');
            
            group = group.setStateFunction('PhaseFlux', WellPhaseFlux(model));
            ctf = ComponentTotalFlux(model);
            ctf.label = 'Q_i';
            group = group.setStateFunction('ComponentTotalFlux', ctf);
            group = group.setStateFunction('ComponentPhaseFlux', WellComponentPhaseFlux(model));
            group = group.setStateFunction('PerforationPressureGradient', PerforationPressureGradient(model));
            group = group.setStateFunction('WellIndex', WellIndex(model));
            group = group.setStateFunction('FacilityWellMapping', FacilityWellMapping(model));
            group = group.setStateFunction('InjectionSurfaceDensity', InjectionSurfaceDensity(model));
            group = group.setStateFunction('ComponentPhaseFractionInjectors', ComponentPhaseFractionInjectors(model));
            group = group.setStateFunction('PerforationMobility', PerforationMobility(model));
            group = group.setStateFunction('PerforationComponentPhaseDensity', PerforationComponentPhaseDensity(model));
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
