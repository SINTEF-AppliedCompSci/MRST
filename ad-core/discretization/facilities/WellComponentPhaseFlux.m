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
            crossFlow = false(numel(wc), nph);
            
            perforationDensity = cell(nph, 1);
            for ph = 1:nph
                flowIntoWell(:, ph) = phaseFlux{ph} > 0;
                crossFlow(:, ph) = flowIntoWell(:, ph) & isInjectorPerforation;
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
                    if ~isempty(componentDensity{ph, c})
                        outflow = ~flowIntoWell(:, ph) & isInjectorPerforation;
                        rhoc = componentDensity{ph, c}(wc);
                        % Compute production fluxes
                        source = rhoc.*q;
                        % Compute injection fluxes
                        if isempty(wcomp{ph})
                            source(outflow) = 0;
                        else
                            source(outflow) = wcomp{ph}.*perforationDensity{ph}(outflow);
                        end
                        % Compute cross-flow fluxes
                        
                        cflux{ph, c} = source;
                    end
                end
            end
        end
    end
end