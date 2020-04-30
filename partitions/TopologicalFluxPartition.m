classdef TopologicalFluxPartition < Partition  
    
    methods
        %-----------------------------------------------------------------%
        function partition = TopologicalFluxPartition(varargin)
            partition = partition@Partition(varargin{:});
            partition.postprocess = false;
        end
        
        %-----------------------------------------------------------------%
        function value = compute(partition, model, state, state0, dt, drivingForces) %#ok
            if isa(model, 'WrapperModel')
                model = model.getReservoirModel();
            end
            flux = state.flux(model.operators.internalConn,:);
            % Get topological permutation
            order = getTopologicalPermutation(model.G, flux, 'W', drivingForces.W, 'padWells', partition.wellPadding > 0);
            % Construct blocks
            nc = max(order);
            nb = min(partition.numBlocks, max(order));
            bz = floor(nc/nb);
            bz = repmat(bz, 1, nb);
            missing = nc - sum(bz);
            bz(1:missing) = bz(1:missing) + 1;
            pos = [0, cumsum(bz(1:end-1))] + 1;
            value = sum(order >= pos,2);
        end
    end
    
end