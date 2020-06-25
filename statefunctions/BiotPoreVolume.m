classdef BiotPoreVolume < StateFunction
    properties
    end
    
    methods
        function gp = BiotPoreVolume(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'BasePoreVolume', 'Dilatation'}, 'BiotPropertyFunctions');
        end
        
        function pv = evaluateOnDomain(prop, model, state)
            vols = model.G.cells.volumes;
            alpha = model.rock.alpha;
            
            [pv, divu] = model.getProps(state, 'BasePoreVolume', 'Dilatation');
            pv = pv + vols.*alpha.*divu;
            
        end
    end
end