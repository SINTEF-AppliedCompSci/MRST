classdef InjectionSurfaceDensity < StateFunction
    % Get injection surface density
    properties

    end
    
    methods
        function gp = InjectionSurfaceDensity(varargin)
            gp@StateFunction(varargin{:});
            gp.dependsOn('FacilityWellMapping');
        end
        function rhoS = evaluateOnDomain(prop, facility, state)
            model = facility.ReservoirModel;
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            W = map.W;
            rhoS = model.getSurfaceDensities();
            % We take the surface density for the first well cell,
            % regardless of active or inactive status for that
            % perforation.
            topcell = arrayfun(@(x) x.cells(1), W);
            reg = model.FlowPropertyFunctions.Density.regions;
            rhoS = rhoS(reg(topcell), :);
            if isfield(W, 'rhoS')
                % Surface density is given on a per-well-basis for the
                % injectors
                rhoS(map.isInjector, :) = vertcat(W(map.isInjector).rhoS);
            end
            nph = size(rhoS, 2);
            tmp = cell(1, nph);
            for i = 1:nph
                tmp{i} = rhoS(:, i);
            end
            rhoS = tmp;
        end
    end
end