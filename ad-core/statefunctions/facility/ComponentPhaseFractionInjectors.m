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

            majorComponent = zeros(1, nph);
            notConcentration = cellfun(@(x) ~x.isConcentration, comp);

            % TODO: the concept of majorComponent can be a little
            % problematic
            for c = 1:ncomp
                if notConcentration(c)
                    majorComponent(comp{c}.phaseIndex) = c;
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
