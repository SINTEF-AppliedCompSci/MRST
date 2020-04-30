classdef UpscaledPartition < Partition

    properties
        parent
        coarseModel
    end
    
    methods
        %-----------------------------------------------------------------%
        function partition = UpscaledPartition(parent, model, partition0)
            partition        = partition@Partition();
            partition.parent = parent;
            coarseModel      = upscaleModelTPFA(model, partition0);
            partition.coarseModel = coarseModel;
        end
        
        %-----------------------------------------------------------------%
        function value = compute(partition, model, state, state0, dt, drivingForces)
            if isa(model, 'WrapperModel')
                model = model.getReservoirModel();
            end
            state  = upscaleState(partition.coarseModel, model, state);
            state0 = upscaleState(partition.coarseModel, model, state0);
            schedule = struct();
            schedule.control = drivingForces;
            schedule = upscaleSchedule(partition.coarseModel, schedule);
            value = partition.parent.compute(partition.coarseModel, state, state0, dt, schedule.control);
            value = value(partition.coarseModel.G.partition);
        end
    end
    
end