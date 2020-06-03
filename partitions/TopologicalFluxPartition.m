classdef TopologicalFluxPartition < Partition  
    % Topological partition based on intercell flux graph
    
    methods
        %-----------------------------------------------------------------%
        function partition = TopologicalFluxPartition(varargin)
            require matlab_bgl
            partition = partition@Partition(varargin{:});
            % Postprocessing may ruin topological order
            partition.postprocess = false;
        end
        
        %-----------------------------------------------------------------%
        function value = compute(partition, model, state, state0, dt, drivingForces) %#ok
            % Compute topological flux partition using matlab_bgl
            if isa(model, 'WrapperModel')
                model = model.getReservoirModel();
            end
            flux = state.flux(model.operators.internalConn,:);
            % Get topological permutation
            order = getTopologicalCellPermutation(model.G, flux, ...
                                                  'W'       , drivingForces.W          , ...
                                                  'padWells', partition.wellPadding > 0);
            % Construct blocks
            nc = max(order);
            nb = min(partition.numBlocks, max(order));
            bz = floor(nc/nb);
            bz = repmat(bz, 1, nb);
            missing = nc - sum(bz);
            bz(1:missing) = bz(1:missing) + 1;
            pos   = [0, cumsum(bz(1:end-1))] + 1;
            value = sum(order >= pos,2);
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
