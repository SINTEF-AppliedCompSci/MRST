classdef MomentumEquations < StateFunction
    
    methods
        function gp = MomentumEquations(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'u', 'lambmech'}, 'state');
            gp = gp.dependsOn('BiotGradP');
            gp.label = 'momentum eqs';
        end
        
        function momentumeqs = evaluateOnDomain(prop, model, state)
            
            G   = model.G;
            op  = model.operators;
            B   = op.B;
            rhs = op.rhs;
            
            u         = model.getProp(state, 'u');
            biotgradp = model.getProp(state, 'BiotGradP');
            lambda    = model.getProp(state, 'lambdamech');
            
            momeqs{1} = B{1,1}*u + B{1,2}*lambda + biotgradp - rhs{1};
            momeqs{2} = B{2,1}*u + B{2,2}*lambda - rhs{2};
    
            momentumeqs = vertcat(momeqs{:});

        end

    end
end

