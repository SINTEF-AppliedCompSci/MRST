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
            isInjectorPerforation = isInj(map.perf2well);
            cflux = cell(nph, ncomp);
            
            % Figure out if things are going into or out from well
            flowIntoWell = false(numel(wc), nph);
            
            perforationDensity = cell(nph, 1);
            for ph = 1:nph
                flowIntoWell(:, ph) = phaseFlux{ph} < 0;
                if any(isInj)
                    perforationDensity{ph} = massDensity{ph}(wc);
                end
            end
            for c = 1:ncomp
                component = model.Components{c};
                if any(isInj)
                    wcomp = component.getPhaseCompositionWell(model, state, W(isInj));
                end
                for ph = 1:nph
                    q = phaseFlux{ph};
                    if ~isempty(componentDensity{c, ph})
                        allowCrossFlow = true;
                        outflow = ~flowIntoWell(:, ph);
                        if allowCrossFlow
                            outflow = outflow & isInjectorPerforation;
                            crossflow = outflow & ~isInjectorPerforation;
                        else
                            crossflow = outflow & false;
                        end
                        rhoc = componentDensity{c, ph}(wc);
                        % Compute production fluxes
                        source = rhoc.*q;
                        % Compute injection fluxes
                        if isempty(wcomp{ph})
                            source(outflow) = 0;
                        else
                            source(outflow) = wcomp{ph}.*perforationDensity{ph}(outflow).*q(outflow);
                        end
                        % Compute cross-flow fluxes
                        if any(crossflow)
                            warning('Crossflow not yet implemented in %s', class(prop));
                        end
                        cflux{c, ph} = source;
                    end
                end
            end
        end
    end
end