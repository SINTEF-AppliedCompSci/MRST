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
            perforationSign = sign(value(phaseFlux));
            crossflow = false(numel(wc), nph);
            totalWellInFlux = zeros(numel(W), nph);
            
            perforationDensity = cell(nph, 1);
            switched_well = false(numel(W), 1);
            % Compute flags for flow into well and cross-flow indicators.
            
            for ph = 1:nph
                crossflow(:, ph) = perforationSign(:, ph) >= 0 & ~isInjectorPerforation;
                all_perf_switched = accumarray(map.perf2well, crossflow(:, ph)) == ncell;
                switched_well = switched_well | all_perf_switched;
            end
            % If all perforations have cross-flow, we are dealing with a
            % switched injector. This routine is responsible for computing
            % component source terms for a given set of phase rates. We
            % assume that the sign of the well has switched since well
            % management is above this routine's paygrade.
            isInj = isInj | switched_well;
            crossflow(switched_well(map.perf2well), :) = false;
            for ph = 1:nph
                if any(isInj)
                    perforationDensity{ph} = massDensity{ph}(wc);
                end
                % If we have cross-flow, we compute total influx for well
                if any(crossflow(:, ph))
                    qi = -min(value(phaseFlux{ph}), 0);
                    totalWellInFlux(:, ph) = accumarray(map.perf2well, qi);
                end
            end
            for c = 1:ncomp
                component = model.Components{c};
                if any(isInj)
                    wcomp = component.getPhaseComponentFractionWell(model, state, W);
                end
                for ph = 1:nph
                    q = phaseFlux{ph};
                    if ~isempty(componentDensity{c, ph})
                        outflow = perforationSign(:, ph) > 0;
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
                            % warning('Crossflow occuring in %d perforations', sum(xflow));
                        end
                        cflux{c, ph} = source;
                    end
                end
            end
        end
    end
end