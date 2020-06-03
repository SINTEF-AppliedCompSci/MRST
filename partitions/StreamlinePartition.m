classdef StreamlinePartition < Partition
    % Partition based on streamlines
    
    properties
        seeds % Streamline seed cells
    end
    
    methods
        %-----------------------------------------------------------------%
        function partition = StreamlinePartition(model, W, varargin)
            assert(any(strcmpi(model.G.type, 'cartGrid')), ...
                'StreamlinePartition requires a Cartesian grid');
            partition = partition@Partition(varargin{:});
            [partition, seedopts] = merge_options(partition, varargin{:});
            if isempty(partition.seeds)
                partition.seeds = partition.getSeedCells(model.G, W, seedopts{:});
            end
            partition.process = true;
            partition.wellPadding = 0.05*(max(model.G.cells.centroids) - min(model.G.cells.centroids));
            partition.wellPadStrategy = 'geometric';
        end
        
        %-----------------------------------------------------------------%
        function seeds = getSeedCells(partition, G, W, varargin) %#ok
            % Get seed cells
            opt = struct('nseeds', 8, 'rWell', []);
            opt = merge_options(opt, varargin{:});
            if isempty(opt.rWell)
                dx        = 0.1*(max(G.nodes.coords) - min(G.nodes.coords));
                opt.rWell = max(dx);
            end
            d = (2*pi/opt.nseeds)/2;
            t = linspace(d, 2*pi - d, opt.nseeds)';
            x = G.cells.centroids;
            seeds = [];
            for i = 1:numel(W)
                if W(i).sign <= 0
                    continue
                end
                xs = opt.rWell.*[cos(t), sin(t)] + x(W(i).cells,:);
                [~, cells] = min(pdist2(x, xs), [], 1);
                seeds = [seeds; cells'];
            end
        end
        
        %-----------------------------------------------------------------%
        function value = compute(partition, model, state, state0, dt, drivingForces) %#ok
            % Compute streamline partition
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

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
