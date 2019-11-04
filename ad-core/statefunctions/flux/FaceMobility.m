classdef FaceMobility < StateFunction & UpwindProperty
    properties (Access = protected)
        upwind_name; % Name of state function where upwind flag comes from
    end
    
    methods
        function fm = FaceMobility(model, upwinding, upwind_name)
            if nargin < 3
                upwind_name = 'PhaseUpwindFlag';
            end
            fm@StateFunction(model);
            fm@UpwindProperty(upwinding)
            fm.upwind_name = upwind_name;
            fm = fm.dependsOn(upwind_name);
            fm = fm.dependsOn('Mobility', 'FlowPropertyFunctions');
        end
        
        function fmob = evaluateOnDomain(prop, model, state)
            flag = prop.getEvaluatedDependencies(state, prop.upwind_name);
            mob = model.getProp(state, 'Mobility');
            nph = numel(mob);
            fmob = cell(1, nph);
            for i = 1:nph
                fmob{i} = prop.faceUpstream(model, state, flag{i}, mob{i});
            end
        end
    end
end