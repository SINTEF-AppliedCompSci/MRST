classdef GravityPermeabilityGradientDG < StateFunction
   
    methods
        function gp = GravityPermeabilityGradientDG(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'Density'}); 
        end
        
        function gRhoKdz = evaluateOnDomain(prop, model, state)
            
            rho = prop.getEvaluatedDependencies(state, 'Density');
            
            g = model.getGravityVector;
            [K, r, c] = permTensor(model.rock, model.G.griddim);
            % FIX for n dimensions.
            assert(model.G.griddim == 2)
            gKdz = [sum(K(:,[1,2]).*g(c([1,2])),2),sum(K(:,[3,4]).*g(c([3,4])),2)];
            gKdz = mat2cell(gKdz, model.G.cells.num, ones(1,model.G.griddim));
            gKdz = SpatialVector(gKdz{:});
            gKdz = gKdz(state.cells,:);
            
            nph = numel(rho);
            gRhoKdz = cell(1,nph);
            for i = 1:nph
                gRhoKdz{i} = rho{i}.*gKdz;
            end
        end
        
    end
        
end
    