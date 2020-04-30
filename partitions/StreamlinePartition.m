classdef StreamlinePartition < Partition
    
    properties
        seeds
        nseeds = 8;
        rWell  = 200;
    end
    
    
    methods
        %-----------------------------------------------------------------%
        function partition = StreamlinePartition(model, W, varargin)
            assert(any(strcmpi(model.G.type, 'cartGrid')), ...
                'StreamlinePartition requires a Cartesian grid');
            partition = partition@Partition(varargin{:});
            partition = merge_options(partition, varargin{:});
            if isempty(partition.seeds)
                partition.seeds = partition.getSeedCells(model.G, W);
            end
            partition.process = true;
            partition.wellPadding = 0.05*(max(model.G.cells.centroids) - min(model.G.cells.centroids));
            partition.wellPadStrategy = 'geometric';
        end
        
        %-----------------------------------------------------------------%
        function seeds = getSeedCells(partition, G, W)
            d = (2*pi/partition.nseeds)/2;
            t = linspace(d, 2*pi - d, partition.nseeds)';
            x = G.cells.centroids;
            seeds = [];
            for i = 1:numel(W)
                if W(i).sign <= 0
                    continue
                end
                xs = partition.rWell.*[cos(t), sin(t)] + x(W(i).cells,:);
                [~, cells] = min(pdist2(x, xs), [], 1);
                seeds = [seeds; cells'];
            end
        end
        
        %-----------------------------------------------------------------%
        function value = compute(partition, model, state, state0, dt, drivingForces) %#ok
            if isa(model, 'WrapperModel')
                model = model.getReservoirModel();
            end
            SF = pollock(model.G, state, partition.seeds);
            SR = pollock(model.G, state, partition.seeds, 'reverse', true);
            S = cellfun(@(sf, sr) vertcat(flipud(sr), sf), SF, SR, 'UniformOutput', false);
            x = model.G.cells.centroids;
            d = nan(model.G.cells.num, numel(S));
            for i = 1:numel(S)
                d(:,i) = min(pdist2(x, S{i}),[],2);
            end
            [~, value] = min(d, [], 2);
        end
        
    end

end