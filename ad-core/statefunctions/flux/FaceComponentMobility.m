classdef FaceComponentMobility < StateFunction & UpwindProperty
    properties (Access = protected)
        upwind_name; % Name of state function where upwind flag comes from
    end
    
    methods
        function fm = FaceComponentMobility(model, upwinding, upwind_name)
            if nargin < 3
                upwind_name = 'PhaseUpwindFlag';
            end
            fm@StateFunction(model);
            fm@UpwindProperty(upwinding)
            fm.upwind_name = upwind_name;
            fm = fm.dependsOn(upwind_name);
            fm = fm.dependsOn('ComponentMobility', 'FlowPropertyFunctions');
        end
        
        function mobf = evaluateOnDomain(prop, model, state)
            flag = prop.getEvaluatedDependencies(state, prop.upwind_name);
            mob = model.getProps(state, 'ComponentMobility');
            [ncomp, nph] = size(mob);
            mobf = cell(ncomp, nph);
            for c = 1:ncomp
                for ph = 1:nph
                    if ~isempty(mob{c, ph})
                        mobf{c, ph} = prop.faceUpstream(model, state, flag{ph}, mob{c, ph});
                    end
                end
            end
        end
    end
end