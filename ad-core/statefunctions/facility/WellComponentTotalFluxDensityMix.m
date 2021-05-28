classdef WellComponentTotalFluxDensityMix < StateFunction
    % Component total flux for wells (with treatment for cross-flow)
    properties

    end
    
    methods
        function gp = WellComponentTotalFluxDensityMix(varargin)
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
            production = ~injection;
            crossflow = (injection & ~isInjector) | ...
                        (production & isInjector);
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
                    surfaceComposition(c, :) = model.ReservoirModel.Components{c}.getPhaseComponentFractionWell(model.ReservoirModel, state, W);
                end
                rem = cellfun(@isempty, surfaceComposition);
                [surfaceComposition{rem}] = deal(zeros(nw, 1));
                % Partial density: Units mass of component per volume of
                % total surface volumetric surface stream
                injectionMass = zeros(nw, ncomp);
                phaseCompi = vertcat(W.compi);
                isRate = ~(strcmpi(targetType, 'resv') | strcmpi(targetType, 'bhp'));
                isRateInjector = isRate & map.isInjector;
                surfaceRates = value(state.FacilityState.surfacePhaseRates);
                surfaceRates(isRateInjector, :) = phaseCompi(isRateInjector, :).*targetValue(isRateInjector);
                surfaceMassRates = surfaceRates.*value(rhoS);
                for ph = 1:nph
                    injectionMass = injectionMass + phaseCompi(:, ph).*surfaceMassRates(:, ph).*[surfaceComposition{:, ph}];
                end
                % Total volumetric flux in connections
                volFlux = model.getProp(state, 'PhaseFlux');
%                 volumeFluxTotal = sum(value(volFlux), 2);
                
                resdens = model.ReservoirModel.getProp(state, 'Density');
                resdens = cellfun(@(x) x(map.cells), resdens, 'unif', false);
                dens = value(resdens);
                
                volFluxValue = value(volFlux);
                injVol = max(volFluxValue, 2);
                surfaceDens = psum*(dens.*injVol);
                surfaceDens = surfaceDens./(psum*injVol);
                
                vi = zeros(nw, ncomp);
                totvol = zeros(nw, nph);
                override = 0;
                for ph = 1:ph
                   wellinj = phaseCompi(:, ph).*surfaceMassRates(:, ph).*[surfaceComposition{:, ph}];
                   override = override - wellinj;
                   wellinj = -max(wellinj, 0); % Into wellbore
                   for c = 1:ncomp

                        mi = phase_flux{c, ph};
                        if isempty(mi)
                            continue
                        end
                        qm = value(min(mi, 0)); % Into wellbore
                        resprod = psum*(qm./dens(:, ph));
                        % Volume of each component in phase (mass in phase
                        % divided by volume)
                        added_volume = resprod + wellinj(:, c)./surfaceDens(:, ph);
                        vi(:, c) = vi(:, c) + added_volume;
                        totvol(:, ph) = totvol(:, ph) + added_volume;
                    end
                end
                ci = bsxfun(@rdivide, vi, sum(vi, 2));
                phasecomp = bsxfun(@rdivide, totvol, sum(totvol, 2));
                bad = any(isnan(ci), 2);
                if any(bad)
                    % Producers that are only injecting - we just need some
                    % value because Newton will hopefully fix this for us.
                    override = abs(override);
                    override = max(override, 1e-12);
                    override = bsxfun(@rdivide, override, sum(override, 2));
                    ci(bad, :) = override(bad, :);
                    phasecomp(bad, :) = override(bad, :); % Just guess...
                end
                vt = zeros(sum(replace), 1);
                for i = 1:nph
                    vt = vt + volFlux{i}(replace);
                end
                
                volFluxTotal = 0;
                for i = 1:nph
                    volFluxTotal = volFluxTotal + volFlux{i};
                end
                massFluxTotal = 0;
                for i = 1:nph
                    massFluxTotal = massFluxTotal + phasecomp(map.perf2well, i).*volFluxTotal.*resdens{i};
                end
                ci_perf = ci(map.perf2well, :);
                for i = 1:ncomp
                    v{i}(replace) = massFluxTotal(replace).*ci_perf(replace, i);
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
