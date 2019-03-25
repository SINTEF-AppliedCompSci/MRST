classdef WellComponentPhaseFlux < GridProperty
    properties

    end
    
    methods
        function gp = WellComponentPhaseFlux(varargin)
            gp@GridProperty(varargin{:});
        end
        function componentPhaseFlux = evaluateOnDomain(prop, facility, state)
            % Preliminaries
            model = facility.ReservoirModel;
            ncomp = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();
            
            % Get fluxes and densities + well map needed
            [map, phaseFlux] = facility.getProps(state, 'FacilityWellMapping', 'PhaseFlux');
            componentDensity = model.getProps(state, 'ComponentPhaseDensity');
            wc = map.cells;
            W = map.W;
            nw = numel(W);
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
            if any(perfIsInjector)
            end
            surfaceComposition = cell(nph, ncomp);
            componentPhaseFlux = cell(nph, ncomp);
            for c = 1:ncomp
                % Store well injector composition
                surfaceComposition(:, c) = model.Components{c}.getPhaseComponentFractionWell(model, state, W);
                for ph = 1:nph
                    % Compute production source terms everywhere. We
                    % overwrite the injection/crossflow terms later on.
                    q = phaseFlux{ph};
                    if ~isempty(componentDensity{c, ph})
                        rhoc = componentDensity{c, ph}(wc);
                        % Compute production fluxes
                        componentPhaseFlux{c, ph} = rhoc.*q;
                    end
                end
            end
            
            
            if any(perfIsInjector)
                massDensity = model.getProp(state, 'Density');
                injPerforationDensity = cellfun(@(x) x(wc(perfInjecting)), massDensity, 'UniformOutput', false);
                
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
                    compi = crossFlowMixture(cflux, compi, map);
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