classdef Partition
   
    properties
        value           = []
        minBlockSize    = -inf
        wellPadding     = 1
        wellPadStrategy = 'topological'
        process         = false
    end
    
    methods
        %-----------------------------------------------------------------%
        function p = Partition(model, varargin)
            p = merge_options(p, varargin{:});
            if p.minBlockSize < 0
                if isa(model, 'WrapperModel')
                    model = model.getReservoirModel();
                end
                p.minBlockSize = model.G.cells.num/20;
            end
        end
        
        %-----------------------------------------------------------------%
        function partition = get(partition, model, state, state0, dt, drivingForces)
            if ~isempty(partition.value) && ~partition.doCompute(state, state0, dt, drivingForces)
                return
            end
            val = partition.compute(model, state, state0, dt, drivingForces);
            if partition.wellPadding > 0
                val = partition.padWells(model, val, drivingForces);
            end
            if partition.process
                val = processPartition(model.G, val);
            end
            partition.value = val;
        end
        
        %-----------------------------------------------------------------%
        function partition = compute(p, model, state, state0, dt, drivingForces) %#ok
            error('Base class not meant for direct use!');
        end
        
        %-----------------------------------------------------------------%
        function ok = doCompute(p, model, state, state0, dt, drivingForces) %#ok
            ok = true;
        end
        
        %-----------------------------------------------------------------%
        function val = padWells(dp, model, val, drivingForces)
            W = drivingForces.W;
            if isempty(W)
                return
            end
            if isa(model, 'WrapperModel')
                model = model.getReservoirModel();
            end
            switch dp.wellPadStrategy
                case 'topological'
                    M = getConnectivityMatrix(model.operators.N);
                    for i = 1:numel(W)
                        cells = false(model.G.cells.num,1);
                        cells(W(i).cells) = true;
                        cell0 = find(cells);
                        for j = 1:dp.wellPadding
                            cells = cells | M*cells;
                        end
                        val(cells) = val(cell0(1));
                    end
                case 'geometric'
                    error('Not implemented yet!')
            end
            val = compressPartition(val);
        end
        
    end
    
end
