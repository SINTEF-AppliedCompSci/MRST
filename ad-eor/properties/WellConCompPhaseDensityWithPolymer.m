classdef WellConCompPhaseDensityWithPolymer < StateFunction
    % Component density to used for each well connection
    properties

    end
    
    methods
        function gp = WellConCompPhaseDensityWithPolymer(varargin)
            gp = gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'ComponentPhaseDensity'}, 'FlowPropertyFunctions');
            % TODO: right the follow?
            gp = gp.dependsOn({'PolymerEffViscMult', 'PolymerViscMult'}, 'PVTPropertyFunctions');
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp.label = '\rho_{wc}';
        end
        
        
        function rhoc = evaluateOnDomain(prop, model, state)
            componentPhaseDensity = model.ReservoirModel.getProps(state, 'ComponentPhaseDensity');
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            function v = ex(x, cells)
                v = [];
                if ~isempty(x) 
                    v=x(cells); 
                end
            end
            f = @(x) ex(x, map.cells);
            rhoc = cellfun(f, componentPhaseDensity, 'UniformOutput', false);
            
            for c = 1 : model.getNumberOfComponents()
                comp = model.ReservoirModel.Components{c};
                if strcmp(comp.name,'polymer')
                    wIx = comp.phaseIndex;
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
