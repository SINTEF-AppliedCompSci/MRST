classdef FacilityFlowDiscretization < StateFunctionGrouping
    properties
        PhaseFlux               % Phase flux in well-bore and to reservoir
        ComponentTotalFlux      % Total mass flux of each component
        ComponentPhaseFlux      % Mass flux of each component, in each phase
        ComponentPhaseDensity   % Density of component in connecting cells
        PressureGradient        % Discrete pressure gradient into well bore
        WellIndex               % Well connection factor / well index
        FacilityWellMapping     % Various mappings used to set up wells
        InjectionSurfaceDensity % Density of injected fluid at surface
        Mobility                % Phase mobilities in connecting cells
    end
    
    methods
        function group = FacilityFlowDiscretization(model)
            group@StateFunctionGrouping('FacilityFluxProps');
            
            group = group.setStateFunction('PhaseFlux', WellPhaseFlux(model));
            ctf = ComponentTotalFlux(model);
            ctf.label = 'Q_i';
            group = group.setStateFunction('ComponentTotalFlux', ctf);
            group = group.setStateFunction('ComponentPhaseFlux', WellComponentPhaseFlux(model));
            group = group.setStateFunction('PressureGradient', PerforationPressureGradient(model));
            group = group.setStateFunction('WellIndex', WellIndex(model));
            group = group.setStateFunction('FacilityWellMapping', FacilityWellMapping(model));
            group = group.setStateFunction('InjectionSurfaceDensity', InjectionSurfaceDensity(model));
            group = group.setStateFunction('Mobility', PerforationMobility(model));
            group = group.setStateFunction('ComponentPhaseDensity', PerforationComponentPhaseDensity(model));
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
