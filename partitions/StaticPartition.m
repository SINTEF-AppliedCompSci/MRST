classdef StaticPartition < Partition
    % Static partition based on precomputed partition vector
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function sp = StaticPartition(partition)
            sp = sp@Partition('value', partition);
        end
        
        %-----------------------------------------------------------------%
        function ok = doCompute(p, model, state, state0, dt, drivingForces) %#ok
            ok = false;
        end
        
    end
    
end