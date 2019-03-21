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
            wellTotalFlux = accumarray(map.perf2well, perfTotalFlux);
            
            wellIsInjector = map.isInjector;
            perfIsInjector = wellIsInjector(map.perf2well);

            perforationSign = sign(perfTotalFlux);
            % Compute flags for flow into well and cross-flow indicators.
            ncell = arrayfun(@(x) numel(x.cells), W);
            perfInjecting = perfTotalFlux > 0;
            perfProducing = perfTotalFlux < 0;
            
            % Compute cross-flow indicators
            crossflow = (perfInjecting &  perfIsInjector) | ...
                        (perfProducing & ~perfIsInjector);

            % If all perforations have cross-flow, we are dealing with a
            % switched injector. This routine is responsible for computing
            % component source terms for a given set of phase rates. We
            % assume that the sign of the well has switched since well
            % management is above this routine's paygrade.
            switched_well = accumarray(map.perf2well, crossflow) == ncell;
            wellIsInjector = wellIsInjector | switched_well;
            crossflow(switched_well(map.perf2well)) = false;
            
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
                wellBoreComposition = cell(1, nph);
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
            
            
            
            if false
                for c = 1:ncomp
                    for ph = 1:nph
                        q = phaseFlux{ph};
                        if ~isempty(componentDensity{c, ph})
                            outflow = perforationSign(:, ph) > 0;
                            rhoc = componentDensity{c, ph}(wc);
                            % Compute production fluxes
                            source = rhoc.*q;
                            % Compute injection fluxes
                            if any(wellIsInjector)
                                if isempty(surfaceComposition{c, ph})
                                    source(outflow) = 0;
                                else
                                    wi = wcomp{ph}(map.perf2well(outflow));
                                    source(outflow) = wi.*perforationDensity{ph}(outflow).*q(outflow);
                                end
                            end
                            % Compute cross-flow fluxes
                            if any(crossflow)
                                % Compute mass inflow in all wells
                                M = accumarray(map.perf2well, -min(value(source), 0));
                                % Mixture density for crossflow in wellbore is
                                % mass flux into well, divided by volume
                                % flowing into the well.
                                rhoMix = M./max(totalWellInFlux(:, ph), 1e-12);
                                source(xflow) = rhoMix(map.perf2well(xflow)).*q(xflow);
                                % warning('Crossflow occuring in %d perforations', sum(xflow));
                            end
                            componentPhaseFlux{c, ph} = source;
                        end
                    end
                end
            end
        end
    end
end