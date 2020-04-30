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