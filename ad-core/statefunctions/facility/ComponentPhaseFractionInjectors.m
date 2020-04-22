classdef ComponentPhaseFractionInjectors < StateFunction
    % Component injection fraction in each phase for wells
    properties

    end

    methods
        function gp = ComponentPhaseFractionInjectors(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn('FacilityWellMapping');
            gp.label = 'Q_{i,\alpha}';
        end
        function compi = evaluateOnDomain(prop, facility, state)
            model = facility.ReservoirModel;
            map = facility.getProp(state, 'FacilityWellMapping');
            ncomp = model.getNumberOfComponents();
            nph = model.getNumberOfPhases();
            surfaceComposition = cell(ncomp, nph);
            W = map.W;
            nw = numel(W);
            for ph = 1:nph
                for c = 1:ncomp
                    % Store well injector composition
                    surfaceComposition(c, :) = model.Components{c}.getPhaseComponentFractionInjection(model, state, W);
                end
            end
            rem = cellfun(@isempty, surfaceComposition);
            [surfaceComposition{rem}] = deal(zeros(nw, 1));
            comp = model.Components;
            isImmiscible = cellfun(@(x) x.isImmiscible, comp);
            mix = any(isImmiscible) && any(~isImmiscible);
            majorComponent = zeros(1, nph);

            notConcentration = cellfun(@(x) ~x.isConcentration, comp);

            if mix
                for c = 1:ncomp
                    % TODO: how to coordinate the isImmiscible and
                    % isConcentration will be a problem here
                    % TODO: only immiscible component can carry other
                    % component? I do not get the full logic here yet
                    if isImmiscible(c)
                        majorComponent(comp{c}.phaseIndex) = c;
                    end
                end
            end
            compi = cell(1, nph);
            for ph = 1:nph
                ci = [surfaceComposition{:, ph}];
                mc = majorComponent(ph);
                if mc > 0
                    otherindex = [1:(mc-1), mc+1:ncomp];
                    ci(:, mc) = ci(:, mc) - sum(ci(:, otherindex(notConcentration(otherindex))), 2);
                end
                compi{ph} = ci;
            end
        end
    end
end
