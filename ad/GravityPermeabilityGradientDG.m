classdef GravityPermeabilityGradientDG < StateFunction
   
    properties
        S
    end
    
    methods
        function gp = GravityPermeabilityGradientDG(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'Density'});
            model = varargin{1};
            dim = model.G.griddim;
            [ii, jj] = blockDiagIndex(dim*ones(dim,1), ones(dim,1));
            gp.S = sparse(ii, jj, 1);
        end
        
        function gRhoKdz = evaluateOnDomain(gp, model, state)
            
            rho = gp.getEvaluatedDependencies(state, 'Density');
            
            g = model.getGravityVector;
            [K, ~, c] = permTensor(model.rock, model.G.griddim);
            gKdz = (K.*g(c))*gp.S;
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
    