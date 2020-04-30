classdef Partition
   
    properties
        value           = []
        numBlocks       = 20;
        wellPadding     = 1
        wellPadStrategy = 'topological'
        process         = false
        postprocess     = true;
        updateFrequency = 1;
    end
    
    methods
        %-----------------------------------------------------------------%
        function partition = Partition(varargin)
            partition = merge_options(partition, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function partition = get(partition, model, state, state0, dt, drivingForces)
            if ~isempty(partition.value) && ~partition.doCompute(state, state0, dt, drivingForces)
                return
            end
            val = partition.compute(model, state, state0, dt, drivingForces);
            if ~partition.postprocess
                partition.value = val;
                return
            end 
            if partition.wellPadding > 0
                val = partition.padWells(model, val, drivingForces);
            end
            if partition.process
                val = processPartition(model.G, val);
            end
            partition.value = val;
        end
        
        %-----------------------------------------------------------------%
        function value = compute(partition, model, state, state0, dt, drivingForces) %#ok
            error('Base class not meant for direct use!');
        end
        
        %-----------------------------------------------------------------%
        function ok = doCompute(partition, model, state, state0, dt, drivingForces) %#ok
            ok = true;
        end
        
        %-----------------------------------------------------------------%
        function val = padWells(partition, model, val, drivingForces)
            W = drivingForces.W;
            if isempty(W)
                return
            end
            if isa(model, 'WrapperModel')
                model = model.getReservoirModel();
            end
            wval = zeros(model.G.cells.num, numel(W));
            np = accumarray(val, 1);
            x = model.G.cells.centroids;
            M = getConnectivityMatrix(model.operators.N);
            for i = 1:numel(W)
                val0 = val(W(i).cells);
                switch partition.wellPadStrategy
                    case 'topological'
                        cells = false(model.G.cells.num,1);
                        cells(W(i).cells) = true;
                        for j = 1:partition.wellPadding
                            cells = cells | M*cells;
                        end
                    case 'geometric'
                        cells = any(pdist2(x, x(W(i).cells,:)) < partition.wellPadding, 2);
                end
                wval(cells, i) = min(val0);
            end
            wval = max(wval,[],2);
            wcells = any(wval > 0, 2);
            val(wcells) = wval(wcells);
            val = compressPartition(val);
        end
    end
    
end
