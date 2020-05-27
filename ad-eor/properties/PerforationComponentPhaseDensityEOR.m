classdef PerforationComponentPhaseDensityEOR < PerforationComponentPhaseDensity
    % Component density to used for each well connection
    properties

    end
    
    methods
        function gp = PerforationComponentPhaseDensityEOR(varargin)
            gp = gp@PerforationComponentPhaseDensity(varargin{:});
            % TODO: right the follow?
            gp = gp.dependsOn({'PolymerEffViscMult', 'PolymerViscMult'}, 'PVTPropertyFunctions');
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp.label = '\rho_{wc}';
        end
        
        
        function rhoc = evaluateOnDomain(prop, model, state)
            rhoc = evaluateOnDomain@PerforationComponentPhaseDensity(prop, model, state);
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            for c = 1 : model.getNumberOfComponents()
                comp = model.ReservoirModel.Components{c};
                if strcmp(comp.name, 'polymer')
                    wIx = 1;
                    % model.ReservoirModel
                    [effviscmult, pviscmult] = model.ReservoirModel.getProps(state, 'PolymerEffViscMult', 'PolymerViscMult');
                    effvismultw = effviscmult(map.cells);
                    pviscmultw = pviscmult(map.cells);
                    rhoc{c, wIx} = rhoc{c, wIx} .* effvismultw ./ pviscmultw;
                end
            end
        end
    end
    
end
