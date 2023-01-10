classdef WellboreComponentPhaseFlux < StateFunction
    
    methods
        %-----------------------------------------------------------------%
        function cf = WellboreComponentPhaseFlux(model)

            cf@StateFunction(model);
            cf = cf.dependsOn({'massFlux'}, 'state');
            cf.label = 'V_{i,\alpha}^w';

        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function v = evaluateOnDomain(prop, model, state)
            
            vt = state.massFlux;
            [mass, flag] = model.getProps(state, 'ComponentPhaseMass', 'PhaseUpwindFlag');
            nc = model.getNumberOfComponents();
            nph   = model.getNumberOfPhases();
            massT = 0;
            for c = 1:nc
                for ph = 1:nph
                    m = mass{c,ph};
                    if isempty(m), continue; end
                    massT = massT + m;
                end
            end
            
            v = cell(nc, nph);
            for c = 1:nc
                for ph = 1:nph
                    m = mass{c,ph};
                    if isempty(m), continue; end
                    x = model.operators.faceUpstr(flag{ph}, m./massT);
                    v{c, ph} = vt.*x;
                end
            end
            
        end
        %-----------------------------------------------------------------%
        
    end
    
end