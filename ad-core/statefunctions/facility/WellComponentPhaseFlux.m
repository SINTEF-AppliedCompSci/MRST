classdef WellComponentPhaseFlux < StateFunction
    % Component flux in each phase for wells
    properties

    end

    methods
        function gp = WellComponentPhaseFlux(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping', 'PhaseFlux', 'ComponentPhaseFractionInjectors'});
            gp = gp.dependsOn({'Density'}, 'PVTPropertyFunctions');
            gp = gp.dependsOn({'WellConCompPhaseDensity'});
            gp.label = 'Q_{i,\alpha}';
        end
        function componentPhaseFlux = evaluateOnDomain(prop, facility, state)
            % Preliminaries
            model = facility.ReservoirModel;
            ncomp = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();

            % Get fluxes and densities + well map needed
            [map, phaseFlux] = prop.getEvaluatedDependencies(state, 'FacilityWellMapping', 'PhaseFlux');
            wellConCompPhaseDensity = model.getProps(state, 'WellConCompPhaseDensity');
            sc = facility.getProps(state, 'ComponentPhaseFractionInjectors');
            wc = map.cells;
            nperf = numel(wc);

            perfTotalFlux = sum(value(phaseFlux), 2);

            wellIsInjector = map.isInjector;
            perfIsInjector = wellIsInjector(map.perf2well);

            % Compute flags for flow into well and cross-flow indicators.
            perfInjecting = perfTotalFlux > 0;

            % If all perforations have cross-flow, we are dealing with a
            % switched injector. This routine is responsible for computing
            % component source terms for a given set of phase rates. We
            % assume that the sign of the well has switched since well
            % management is above this routine's paygrade.


            % Get phase density if we are injecting
            componentPhaseFlux = cell(ncomp, nph);
            for c = 1:ncomp
                % Store well injector composition
                for ph = 1:nph
                    % Compute production source terms everywhere. We
                    % overwrite the injection/crossflow terms later on.
                    q = phaseFlux{ph};
                    if ~isempty(wellConCompPhaseDensity{c, ph})
                        rhoc = wellConCompPhaseDensity{c, ph};
                        componentPhaseFlux{c, ph} = rhoc.*q;
                    end
                end
            end
            if any(perfIsInjector)
                massDensity = model.getProp(state, 'Density');
                injPerforationDensity = cellfun(@(x) x(wc(perfInjecting)), massDensity, 'UniformOutput', false);

                for ph = 1:nph
                    q = phaseFlux{ph};

                    cflux = zeros(nperf, ncomp);
                    for c = 1:ncomp
                        v = componentPhaseFlux{c, ph};
                        if ~isempty(v)
                            cflux(:, c) = value(v);
                        end
                    end

                    compi = sc{ph};
                    comp = model.Components;
                    notconc = cellfun(@(c)(~c.isConcentration), comp);
                    compi(:, notconc) = compi(:, notconc)./max(sum(compi(:, notconc),2), 1e-10);
                    compi(:, notconc) = crossFlowMixture(cflux(:, notconc), compi(:, notconc), map);
                    for c = 1:ncomp
                        if ~isempty(componentPhaseFlux{c, ph})
                            perfCompDens = compi(map.perf2well(perfInjecting), c).*injPerforationDensity{ph};
                            componentPhaseFlux{c, ph}(perfInjecting) = q(perfInjecting).*perfCompDens;
                        end
                    end
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
