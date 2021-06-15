classdef WellComponentTotalVolumeBalanceCrossflow < StateFunction
    % Component total flux for wells (with treatment for cross-flow)
    properties
        onlyLocalDerivatives = true;
    end
    
    methods
        function gp = WellComponentTotalVolumeBalanceCrossflow(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp = gp.dependsOn({'ComponentPhaseFlux', 'InjectionSurfaceDensity'});
        end
        
        function v = evaluateOnDomain(prop, model, state)
            ncomp = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();
            v = cell(ncomp, 1);
            phase_flux = prop.getEvaluatedDependencies(state, 'ComponentPhaseFlux');
            
            for c = 1:ncomp
                % Loop over phases where the component may be present
                for ph = 1:nph
                    % Check if present
                    m = phase_flux{c, ph};
                    if ~isempty(m)
                        if isempty(v{c})
                            v{c} = m;
                        else
                            v{c} = v{c} + m;
                        end
                    end
                end
            end
            % If we have cross-flow and/or we are injecting more than one
            % component, we need to ensure that injecting perforation
            % composition mixtures equal the mix entering the wellbore
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            W   = map.W;
            massFlux  = value(v');
            massFluxTotal = sum(massFlux,2);
            psum = map.perforationSum;

            isInjector = map.isInjector(map.perf2well);
            injection  = massFluxTotal > 0;
            production = massFluxTotal < 0;
            crossflow = (injection & ~isInjector) | ...
                        (production & isInjector);% & abs(value(massFluxTotal)) > 1e-8;
            wellHasXFlow = psum*crossflow;
            replace = injection & wellHasXFlow(map.perf2well);

            if any(replace)
                rhoS = prop.getEvaluatedDependencies(state, 'InjectionSurfaceDensity');

                ws = state.wellSol(map.active);
                targetType = reshape(arrayfun(@(x) x.type, ws, 'UniformOutput', false), [], 1);
                targetValue = vertcat(ws.val);
%                 isZero = ~strcmpi(targets, 'bhp') & val == 0;
                nw = numel(W);
                surfaceComposition = cell(ncomp, nph);
                for c = 1:ncomp
                    % Store well injector composition
                    surfaceComposition(c, :) = model.ReservoirModel.Components{c}.getPhaseComponentFractionInjection(model.ReservoirModel, state, W);
                end
                rem = cellfun(@isempty, surfaceComposition);
                [surfaceComposition{rem}] = deal(zeros(nw, 1));
                % Partial density: Units mass of component per volume of
                % total surface volumetric surface stream
                injectionMass = cell(ncomp, nph);
                [injectionMass{:}] = deal(0);
                phaseCompi = vertcat(W.compi);
                isRate = ~(strcmpi(targetType, 'resv') | strcmpi(targetType, 'bhp'));
                isRateInjector = isRate & map.isInjector;
                qs = state.FacilityState.surfacePhaseRates;

                target = vertcat(ws.val);
                for ph = 1:nph
                    q_surf = (~isRateInjector.*value(qs{ph}) + isRateInjector.*target);
%                     q_surf = qs{ph};
                    qs_ph = rhoS{ph}.*phaseCompi(:, ph).*q_surf;
                    for c = 1:ncomp
                        if any(surfaceComposition{c, ph})
                            injectionMass{c, ph} = qs_ph.*surfaceComposition{c, ph};
                        end
                    end
                end
                
                % Total volumetric flux in connections
                volFlux = model.getProp(state, 'PhaseFlux');
                
                resdens = model.ReservoirModel.getProp(state, 'Density');
                resdens = cellfun(@(x) x(map.cells), resdens, 'unif', false);
                
                wellComponentTotalMass = cell(1, ncomp);
                wellComponentTotalMassXflow = cell(1, ncomp);
                wellPhaseTotalMass = cell(1, nph);
                [wellPhaseTotalMass{:}, wellComponentTotalMass{:}, wellComponentTotalMassXflow{:}] = deal(0);
                phaseMassComposition = cell(ncomp, nph);
                for ph = 1:ph
                   for c = 1:ncomp
                        mi = phase_flux{c, ph};
                        if isempty(mi)
                            continue
                        end
                        from_surface = -max(injectionMass{c, ph}, 0); % 
                        qm = min(mi, 0);
                        qm = value(qm);
                        qm = prop.reduce(qm); % Into wellbore
%                         added_mass = psum*qm + from_surface;
                        xflow = psum*qm;
                        added_mass = xflow + from_surface;
                        wellPhaseTotalMass{ph} = wellPhaseTotalMass{ph} + added_mass;
                        wellComponentTotalMassXflow{c} = wellComponentTotalMassXflow{c} + xflow;
                        wellComponentTotalMass{c} = wellComponentTotalMass{c} + added_mass;
                        phaseMassComposition{c, ph} = added_mass;
                   end
                   tot = 0;
                   for c = 1:ncomp
                       x = phaseMassComposition{c, ph};
                       if ~isempty(x)
                           x = max(abs(x), 1e-12);
                           tot = tot + x;
                           phaseMassComposition{c, ph} = x;
                       end
                   end
                   for c = 1:ncomp
                       x = phaseMassComposition{c, ph};
                       if ~isempty(x)
                           phaseMassComposition{c, ph} = x./tot;
                       end
                   end
                end
                
                volFluxTotal = 0;
                for i = 1:nph
                    volFluxTotal = volFluxTotal + volFlux{i};
                end
                phaseVolumeIntoWellbore = cell(1, nph);
                active = map.perf2well(replace);
                avgVolume = false;
                
                tot = 0;
                for ph = 1:nph
                    if avgVolume
                        qf = value(volFlux{ph});
                        r = psum*value(resdens{ph}.*qf)./(psum*qf);
                        r(isnan(r)) = 1;
                        r = r(active);
                    else
                        r = resdens{ph}(replace);
                    end
                    x = wellPhaseTotalMass{ph}(active)./r;
                    x = max(abs(x), 1e-12);
                    tot = tot + x;
                    phaseVolumeIntoWellbore{ph} = x;
                end
                phaseVolumeFractions = cell(1, nph);
                for ph = 1:nph
                    phaseVolumeFractions{ph} = phaseVolumeIntoWellbore{ph}./tot;
                end
                for c = 1:ncomp
                    v{c}(replace) = 0;
                end
                replaceMassFlux = cell(1, nph);
                vft = volFluxTotal(replace);
                for ph = 1:nph
                    rho = resdens{ph}(replace);
                    replaceMassFlux{ph} = phaseVolumeFractions{ph}.*vft.*rho;
                end
                v_r = prop.computeCrossFlux(model, state, 'specialsauce', map, replace, ncomp, nph, resdens, phase_flux, replaceMassFlux, wellComponentTotalMass, wellComponentTotalMassXflow, phaseMassComposition, volFluxTotal, rhoS, phaseCompi, surfaceComposition, qs, injectionMass);
                v_f = prop.computeCrossFlux(model, state, 'phasemassfractions', map, replace, ncomp, nph, resdens, phase_flux, replaceMassFlux, wellComponentTotalMass, wellComponentTotalMassXflow, phaseMassComposition, volFluxTotal, rhoS, phaseCompi, surfaceComposition, qs, injectionMass);
                for c = 1:ncomp
                    ii = map.isInjector(active);
                    new = v_f{c}.*~ii;
                    if any(ii)
                        new = new + v_r{c}.*ii;
                    end
                    v{c}(replace) = new;
                end
            end
        end
        
        function v = computeCrossFlux(prop, model, state, masstype, map, replace, ncomp, nph, resdens, phase_flux, replaceMassFlux, wellComponentTotalMass, ...
                wellComponentTotalMassXflow, phaseMassComposition, volFluxTotal, rhoS, phaseCompi, surfaceComposition, qs, injectionMass)
            psum = map.perforationSum;
            active = map.perf2well(replace);
            vft = volFluxTotal(replace);
            injVol = prop.reduce(volFluxTotal);
            injVol = max(injVol, 0);
            injVolTot = psum*injVol;

            v = cell(1, ncomp);
            switch lower(masstype)
                case 'exactmassfractions'
                    replaceMassFluxTotal = 0;
                    for ph = 1:nph
                        replaceMassFluxTotal = replaceMassFluxTotal + replaceMassFlux{ph};
                    end
                    tot = 0;
                    for c = 1:ncomp
                        tot = tot + wellComponentTotalMass{c};
                    end
                    xi = abs(value(wellComponentTotalMass)) + 1e-12;
                    xi = bsxfun(@rdivide, xi, sum(xi, 2));
                    for c = 1:ncomp
                        x = xi(active, c);
                        v{c} = x.*replaceMassFluxTotal;
                    end
                case 'phasemassfractions'
                    for c = 1:ncomp
                        qi = 0;
                        for ph = 1:nph
                            x = phaseMassComposition{c, ph};
                            if ~isempty(x)
                                flux = replaceMassFlux{ph};
                                qi = qi + flux.*x(active);
                            end
                        end
                        v{c} = qi;
                    end
                case 'exactmass'
                    for c = 1:ncomp
                        mass = -wellComponentTotalMass{c};
                        rho = mass./injVolTot;

                        rp = rho(active);
                        rp = value(rp); % Singular systems if not value?
                        qi = rp.*vft;
                        v{c} = qi;
                    end
                case 'specialsauce'
                    % Next do the well part - assume it matches net
                    % volume injected
                    tmp = prop.reduce(volFluxTotal);
                    if prop.onlyLocalDerivatives
                        delta = -value(volFluxTotal) + volFluxTotal;
                        phases = zeros(numel(map.W), nph, ncomp);
                        for ph = 1:nph
                            for c = 1:ncomp
                                ii = value(phaseMassComposition{c, ph});
                                if ~isempty(ii)
                                    phases(:, ph, c) = ii;
                                end
                            end
                        end
                    else
                        delta = 0;
                    end
                    
                    netVol = psum*tmp;
                    netVol = max(netVol, 0);
                    netVolPerf = netVol(map.perf2well) + delta;
                    injVolTotPerf = injVolTot(map.perf2well) + delta;
                    ratio = netVolPerf./injVolTotPerf;
                    ratio(~isfinite(value(ratio))) = 0;
                    
                    tmp = psum*injVol;
                    w = volFluxTotal./tmp(map.perf2well);
                    w(~isfinite(value(w))) = 1;
                    for c = 1:ncomp
                        % Exact mass-balance for crossflow
                        mass = -wellComponentTotalMassXflow{c}(active);
                        rho = mass./injVolTotPerf(replace);

%                         if prop.onlyLocalDerivatives
%                             rho_mix = 0;
%                             for ph = 1:nph
%                                 rho_mix = rho_mix + phases(map.perf2well, ph, c).*resdens{ph};
%                             end
%                             rho_mix = rho_mix(replace);
%                             rho = rho - value(rho_mix) + rho_mix;
%                         end
                        qi = rho.*vft;
                        % Volume-based approach for remainder
                        extra = 0;
                        qt = volFluxTotal.*ratio;
                        rhoS = model.ReservoirModel.getSurfaceDensities();
%                         rhoS = rhoS(1, :);
                        added = false;
                        for ph = 1:nph
%                             ci = surfaceComposition{c, ph};
%                             tmp = phaseCompi(:, ph).*ci;
%                             add = qs{ph}.*rhoS(ph).*tmp;
                            add = injectionMass{c, ph};
                            if isempty(add) || all(value(add) == 0)
                                continue
                            end
                            add = max(add, 0);
%                             add = value(add);
                            extra = extra + add(map.perf2well).*w;
                            added = true;
%                             extra = extra + tmp(map.perf2well).*resdens{ph}.*qt;
                        end
                        if added
                            v{c} = qi + extra(replace);
                        else
                            v{c} = qi;
                        end
                    end
                otherwise
                    error('%s is not supported', masstype);
            end
        end
        
        function x = reduce(prop, x)
            if prop.onlyLocalDerivatives
                x = value(x);
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
