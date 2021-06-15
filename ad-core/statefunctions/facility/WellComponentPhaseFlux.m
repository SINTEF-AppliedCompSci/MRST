classdef WellComponentPhaseFlux < StateFunction
    % Component flux in each phase for wells
    properties

    end
    
    methods
        function gp = WellComponentPhaseFlux(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping', 'PhaseFlux', 'ComponentPhaseDensity'});
            gp = gp.dependsOn({'Density'}, 'PVTPropertyFunctions');
            gp.label = 'Q_{i,\alpha}';
        end
        function componentPhaseFlux = evaluateOnDomain(prop, facility, state)
            % Preliminaries
            model = facility.ReservoirModel;
            ncomp = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();
            
            % Get fluxes and densities + well map needed
            [map, phaseFlux, componentDensity] = prop.getEvaluatedDependencies(state, 'FacilityWellMapping', 'PhaseFlux', 'ComponentPhaseDensity');
            wc = map.cells;
            W = map.W;
            nw = numel(W);
            componentPhaseFlux = cell(ncomp, nph);
            if nw == 0
                return;
            end
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
            surfaceComposition = cell(ncomp, nph);
            for c = 1:ncomp
                % Store well injector composition
                surfaceComposition(c, :) = model.Components{c}.getPhaseComponentFractionInjection(model, state, W);
                for ph = 1:nph
                    % Compute production source terms everywhere. We
                    % overwrite the injection/crossflow terms later on.
                    q = phaseFlux{ph};
                    rhoc = componentDensity{c, ph};
                    if ~isempty(rhoc)
                        % Compute production fluxes
                        componentPhaseFlux{c, ph} = rhoc.*q;
                    end
                end
            end
            if any(perfIsInjector)
                massDensity = model.PVTPropertyFunctions.get(model, state, 'Density');
                injPerforationDensity = cellfun(@(x) x(wc(perfInjecting)), massDensity, 'UniformOutput', false);
                isMass = ~cellfun(@(c) isa(c, 'ConcentrationComponent'), model.Components);
                rem = cellfun(@isempty, surfaceComposition);
                [surfaceComposition{rem}] = deal(zeros(nw, 1));
                for ph = 1:nph
                    q = phaseFlux{ph};
                    
                    compi = [surfaceComposition{:, ph}];
                    cflux = zeros(nperf, ncomp);
                    for c = 1:ncomp
                        v = componentPhaseFlux{c, ph};
                        if ~isempty(v)
                            cflux(:, c) = value(v);
                        end
                    end
                    compi(:, isMass) = compi(:, isMass)./max(sum(compi(:, isMass),2), 1e-10);
                    compi(:, isMass) = crossFlowMixture(cflux(:, isMass), compi(:, isMass), map);
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
