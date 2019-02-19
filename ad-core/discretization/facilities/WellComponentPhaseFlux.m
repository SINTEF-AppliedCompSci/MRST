classdef WellComponentPhaseFlux < GridProperty
    properties

    end
    
    methods

        function cflux = evaluateOnDomain(prop, facility, state)
            model = facility.ReservoirModel;
            ncomp = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();
            
            map = facility.getProp(state, 'FacilityWellMapping');
            phaseFlux = facility.getProp(state, 'PhaseFlux');
            [componentDensity, massDensity] = model.getProps(state, 'ComponentPhaseDensity', 'Density');
            
            
            isInj = map.isInjector;
            wc = map.cells;
            W = map.W;
            
            ncell = arrayfun(@(x) numel(x.cells), W);
            isInjectorPerforation = isInj(map.perf2well);
            cflux = cell(nph, ncomp);
            
            % Figure out if things are going into or out from well
            flowIntoWell = false(numel(wc), nph);
            crossflow = false(numel(wc), nph);
            totalWellInFlux = zeros(numel(W), nph);
            
            perforationDensity = cell(nph, 1);
            switched_well = false(numel(W), 1);
            
            for ph = 1:nph
                flowIntoWell(:, ph) = phaseFlux{ph} < 0;
                crossflow(:, ph) = ~flowIntoWell(:, ph) & ~isInjectorPerforation;
                all_perf_switched = accumarray(map.perf2well, crossflow(:, ph)) == ncell;
                switched_well = switched_well | all_perf_switched;
            end
            isInj = isInj | switched_well;
            for ph = 1:nph
                if any(isInj)
                    perforationDensity{ph} = massDensity{ph}(wc);
                end
                
                if any(crossflow(:, ph))
                    qi = -min(value(phaseFlux{ph}), 0);
                    totalWellInFlux(:, ph) = accumarray(map.perf2well, qi);
                end
            end
            for c = 1:ncomp
                component = model.Components{c};
                if any(isInj)
                    wcomp = component.getPhaseCompositionWell(model, state, W);
                end
                for ph = 1:nph
                    q = phaseFlux{ph};
                    if ~isempty(componentDensity{c, ph})
                        outflow = ~flowIntoWell(:, ph);
                        rhoc = componentDensity{c, ph}(wc);
                        % Compute production fluxes
                        source = rhoc.*q;
                        % Compute injection fluxes
                        if any(isInj)
                            if isempty(wcomp{ph})
                                source(outflow) = 0;
                            else
                                wi = wcomp{ph}(map.perf2well(outflow));
                                source(outflow) = wi.*perforationDensity{ph}(outflow).*q(outflow);
                            end
                        end
                        % Compute cross-flow fluxes
                        xflow = crossflow(:, ph);
                        if any(xflow)
                            % Compute mass inflow in all wells
                            M = accumarray(map.perf2well, -min(value(source), 0));
                            % Mixture density for crossflow in wellbore is
                            % mass flux into well, divided by volume
                            % flowing into the well.
                            rhoMix = M./max(totalWellInFlux(:, ph), 1e-12);
                            source(xflow) = rhoMix(map.perf2well(xflow)).*q(xflow);
                            warning('Crossflow occuring in %d perforations', sum(xflow));
                        end
                        cflux{c, ph} = source;
                    end
                end
            end
        end
    end
end