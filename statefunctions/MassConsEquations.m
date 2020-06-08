classdef MassConsEquations < StateFunction
    
    methods
        function gp = MomentumEquations(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'pressure', 'lambda'}, 'state');
            gp.label = 'mass cons eqs';
        end
        
        function masseqs = evaluateOnDomain(prop, model, state)
            
            G  = model.G;
            op = model.operators;
            
            B = op.B;
            rhs = op.rhs;
            
            p      = model.getProp(state, 'pressure');
            lambda = model.getProp(state, 'lambda');
            pv     = model.getProp(state, 'PoreVolume');
            
            masseqs{1} = B{1,1}*p + B{1,2}*lambda - rhs{1};
            masseqs{2} = B{2,1}*p + B{2,2}*lambda - rhs{2};
    
            masseqs = vertcat(masseqs{:});

        end

    end
end

