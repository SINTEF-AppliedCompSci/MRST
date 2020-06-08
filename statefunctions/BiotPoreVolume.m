classdef BiotPoreVolume < StateFunction
    properties
    end
    
    methods
        function gp = BiotPoreVolume(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'BasePoreVolume', 'Dilatation'});
        end
        
        function pv = evaluateOnDomain(prop, model, state)
            vols = model.G.cells.volumes;
            alpha = model.rock.alpha;
            
            [pv, divu] = prop.getEvaluatedDependencies(state, 'BasePoreVolume', 'Dilatation');
            pv = pv + vols.*alpha.*divu;
            
        end
    end
end