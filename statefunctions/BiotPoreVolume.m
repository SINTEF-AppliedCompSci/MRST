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
            
            [pv, divu] = model.getProps(state, 'BasePoreVolume', 'Dilatation');
            pv = pv + vols.*divu;
            
        end
    end
end