classdef UpscaledPartition < Partition
    % Upscaled partition wrapped around existing partition

    properties
        parent      % Fine-scale parent model
        coarseModel % Coarse model
    end
    
    methods
        %-----------------------------------------------------------------%
        function partition = UpscaledPartition(parent, model, partition0)
            % Construc upscaled partition
            partition        = partition@Partition();
            partition.parent = parent;
            % Construct upscaled coarse model based on partition0
            partition.coarseModel = upscaleModelTPFA(model, partition0);
        end
        
        %-----------------------------------------------------------------%
        function value = compute(partition, model, state, state0, dt, drivingForces)
            % Compute partition
            if isa(model, 'WrapperModel')
                model = model.getReservoirModel();
            end
            % Upscale state and forces
            state  = upscaleState(partition.coarseModel, model, state);
            state0 = upscaleState(partition.coarseModel, model, state0);
            schedule = struct();
            schedule.control = drivingForces;
            schedule = upscaleSchedule(partition.coarseModel, schedule);
            % Call parent partition to compute
            value = partition.parent.compute(partition.coarseModel, state, state0, dt, schedule.control);
            % Map fo fine-scale model
            value = value(partition.coarseModel.G.partition);
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
