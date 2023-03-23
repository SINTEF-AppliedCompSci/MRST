classdef PoroelasticFractureTransMult < StateFunction
    % Transmissibility mulitplier arising from non-linear permeability
    properties
    end
    
    methods
        function gp = PoroelasticFractureTransMult(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'PoroelasticFracturePerm'});
        end
        function transMult = evaluateOnDomain(prop, model, state)
            N = model.G.faces.neighbors;
            intInx = all(N ~= 0, 2);
            N = N(intInx, :);
            perm_f0 = model.rock.perm;
            perm_f = prop.getEvaluatedDependencies(state, 'PoroelasticFracturePerm');
            transMult = zeros(size(N,1),1);
            harmonic = @(a,b) 1 ./ (1./a + 1./b);
            for n = 1:size(N,1)
                % K_new / K_ref at face
                transMult(n) = harmonic(perm_f(N(n,1)).val, perm_f(N(n,2)).val)... 
                               ./ harmonic(perm_f0(N(n,1)), perm_f0(N(n,2))); 
            end
        end
    end
end