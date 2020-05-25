classdef PerforationMobility < StateFunction
    % Mobility in each perforated cell of a well
    properties
    end
    
    methods
        function gp = PerforationMobility(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping'});
            gp = gp.dependsOn({'Mobility'}, 'FlowPropertyFunctions');
            gp.label = '\lambda_{wc}'; 
        end
        
        function mobw = evaluateOnDomain(prop, model, state)
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            mob = model.ReservoirModel.getProps(state, 'Mobility');
            mobw = cellfun(@(x) x(map.cells), mob, 'UniformOutput', false);            
        end
    end   
end