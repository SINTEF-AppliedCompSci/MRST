classdef Partition
    % Class for generating partitions of a grid based on one or more of the
    % following: model, state, state0, dt, drivingForces
    
    properties
        value           = []            % Partition value, one integer for each cell
        numBlocks       = 20;           % Number of blocks in partition
        postprocess     = true;         % Postprocess partition (pad wells and call processPartition)
        process         = false         % Call processParttition after computing
        wellPadding     = 1             % Padding around wells
        wellPadStrategy = 'topological' % Well padding strategy
    end
    
    methods
        %-----------------------------------------------------------------%
        function partition = Partition(varargin)
            % Constructor
            partition = merge_options(partition, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function partition = get(partition, model, state, state0, dt, drivingForces)
            % Get partition. This function shold not need to be overloaded
            if ~isempty(partition.value) && ~partition.doCompute(state, state0, dt, drivingForces)
                % Partition already computed, and should not be recomputed
                return
            end
            % Compute partition
            val = partition.compute(model, state, state0, dt, drivingForces);
            % Postprocess partition
            if ~partition.postprocess
                partition.value = val;
                return
            end
            % Pad wells
            if partition.wellPadding > 0
                val = partition.padWells(model, val, drivingForces);
            end
            % Process partition
            if partition.process
                val = processPartition(model.G, val);
            end
            partition.value = val;
        end
        
        %-----------------------------------------------------------------%
        function value = compute(partition, model, state, state0, dt, drivingForces) %#ok
            % Main function. Subclasses should overload this function
            error('Base class not meant for direct use!');
        end
        
        %-----------------------------------------------------------------%
        function ok = doCompute(partition, model, state, state0, dt, drivingForces) %#ok
            % Decide if partition should be recomputed
            ok = true;
        end
        
        %-----------------------------------------------------------------%
        function val = padWells(partition, model, val, drivingForces)
            % Pad wells
            W = drivingForces.W;
            if isempty(W)
                return
            end
            if isa(model, 'WrapperModel')
                model = model.getReservoirModel();
            end
            wval = zeros(model.G.cells.num, numel(W));
            x    = model.G.cells.centroids;
            C    = getConnectivityMatrix(model.operators.N);
            for i = 1:numel(W)
                % For each well, add pad of cells
                val0 = val(W(i).cells);
                switch partition.wellPadStrategy
                    case 'topological'
                        cells = false(model.G.cells.num,1);
                        cells(W(i).cells) = true;
                        for j = 1:partition.wellPadding
                            cells = cells | C*cells;
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
